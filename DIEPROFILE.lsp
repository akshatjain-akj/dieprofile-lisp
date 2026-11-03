;;; ===================================================================
;;;  DIEPROFILE.lsp  -  Wire-drawing die profile generator for AutoCAD
;;;
;;;  Command : DIEPROFILE
;;;
;;;  You give : a name (any), the FINAL shaped profile (closed polyline /
;;;             spline), and ALL die areas from INITIAL (round) to FINAL.
;;;  It makes : 1. the round die circle (diameter from the first area)
;;;             2. every intermediate die profile with EXACTLY its area
;;;             3. (optional) a separate 3D solid of the whole process,
;;;                each die section spaced along Z using the die angle.
;;;
;;;  Method   : every point of the final shape moves radially toward the
;;;             round die by a fraction t; corners are rounded a little
;;;             more at each stage; t is solved so the drawn polyline area
;;;             equals the target (error < 0.000001 mm2).  AutoCAD's own
;;;             area of every drawn polyline is read back and printed.
;;; ===================================================================

(vl-load-com)

(setq *DP-NPTS* 720)        ; points per profile (more = smoother)
(setq *DP-TOL*  0.5)        ; % : allowed difference shape area vs final area

;;; ------------------------------------------------------------------
;;;  Geometry helpers  (points are 2D lists (x y), profile is closed)
;;; ------------------------------------------------------------------

;; signed area (shoelace), + = counter-clockwise
(defun dp:area (pts / a lst p q)
  (setq lst (append pts (list (car pts))) a 0.0)
  (while (cdr lst)
    (setq p (car lst) q (cadr lst)
          a (+ a (- (* (car p) (cadr q)) (* (car q) (cadr p))))
          lst (cdr lst)))
  (* 0.5 a))

;; area centroid
(defun dp:centroid (pts / a cx cy lst p q c)
  (setq lst (append pts (list (car pts))) a 0.0 cx 0.0 cy 0.0)
  (while (cdr lst)
    (setq p (car lst) q (cadr lst)
          c (- (* (car p) (cadr q)) (* (car q) (cadr p)))
          a (+ a c)
          cx (+ cx (* (+ (car p) (car q)) c))
          cy (+ cy (* (+ (cadr p) (cadr q)) c))
          lst (cdr lst)))
  (list (/ cx (* 3.0 a)) (/ cy (* 3.0 a))))

;; n points equally spaced along a closed curve
(defun dp:sample (ent n / len k p pts)
  (setq len (vlax-curve-getDistAtParam ent (vlax-curve-getEndParam ent))
        k 0)
  (repeat n
    (setq p   (vlax-curve-getPointAtDist ent (* len (/ (float k) n)))
          pts (cons (list (car p) (cadr p)) pts)
          k   (1+ k)))
  (reverse pts))

(defun dp:firstn (l k / r)
  (repeat k (setq r (cons (car l) r) l (cdr l)))
  (reverse r))

(defun dp:lastn (l k) (reverse (dp:firstn (reverse l) k)))

(defun dp:normang (d)
  (while (> d pi) (setq d (- d (* 2.0 pi))))
  (while (<= d (- pi)) (setq d (+ d (* 2.0 pi))))
  d)

;; Matching point on the round die for every shape point (ray from the
;; shape centroid g).  Angles are forced to increase steadily so the
;; profiles never fold, even for Z-shapes with deep notches.
(defun dp:circpts (pts g cc r / prev a un m fw bw)
  (foreach p pts
    (setq a (angle g p))
    (if prev (setq a (+ prev (dp:normang (- a prev)))))
    (setq un (cons a un) prev a))
  (setq un (reverse un) m -1e99)
  (foreach a un (setq m (max m a) fw (cons m fw)))
  (setq fw (reverse fw) m 1e99)
  (foreach a (reverse un) (setq m (min m a) bw (cons m bw)))
  (mapcar '(lambda (f b) (polar cc (/ (+ f b) 2.0) r)) fw bw))

;; closed moving-average over 2m+1 points
(defun dp:box (pts m n / w ext lead tail sx sy res)
  (setq w    (+ m m 1)
        ext  (append (dp:lastn pts m) pts (dp:firstn pts m))
        lead ext tail ext sx 0.0 sy 0.0)
  (repeat w
    (setq sx (+ sx (caar lead)) sy (+ sy (cadar lead)) lead (cdr lead)))
  (repeat n
    (setq res (cons (list (/ sx w) (/ sy w)) res))
    (if lead
      (setq sx (+ sx (- (caar lead) (caar tail)))
            sy (+ sy (- (cadar lead) (cadar tail)))
            lead (cdr lead)))
    (setq tail (cdr tail)))
  (reverse res))

;; corner rounding (3 passes of moving average ~ gaussian)
(defun dp:smooth (pts m n)
  (if (>= m 1) (repeat 3 (setq pts (dp:box pts m n))))
  pts)

;; point-by-point blend  shape -> circle
(defun dp:blend (ps cs tt)
  (mapcar '(lambda (p c)
             (list (+ (car p)  (* tt (- (car c)  (car p))))
                   (+ (cadr p) (* tt (- (cadr c) (cadr p))))))
          ps cs))

;; find t in [0,1] so that area(blend) = target  (bisection)
(defun dp:solve (ps cs target / lo hi mid)
  (setq lo 0.0 hi 1.0)
  (repeat 50
    (setq mid (/ (+ lo hi) 2.0))
    (if (< (dp:area (dp:blend ps cs mid)) target)
      (setq lo mid)
      (setq hi mid)))
  (/ (+ lo hi) 2.0))

;; full profile for one target area -> list of points
(defun dp:profile (P C target n rf / t0 m ps)
  (setq t0 (dp:solve P C target)
        m  (fix (+ 0.5 (* t0 n rf))))
  (if (> m (/ n 4)) (setq m (/ n 4)))
  (setq ps (dp:smooth P m n))
  (dp:blend ps C (dp:solve ps C target)))

;;; ------------------------------------------------------------------
;;;  Drawing helpers
;;; ------------------------------------------------------------------
(defun dp:pline (pts lay z)
  (entmakex
    (append
      (list '(0 . "LWPOLYLINE") '(100 . "AcDbEntity") (cons 8 lay)
            '(100 . "AcDbPolyline") (cons 90 (length pts)) '(70 . 1)
            (cons 38 z))
      (mapcar '(lambda (p) (cons 10 p)) pts))))

(defun dp:circle (cc r lay z)
  (entmakex (list '(0 . "CIRCLE") (cons 8 lay)
                  (list 10 (car cc) (cadr cc) z) (cons 40 r))))

(defun dp:text (pt h str lay)
  (entmakex (list '(0 . "TEXT") (cons 8 lay)
                  (list 10 (car pt) (cadr pt) 0.0)
                  (cons 40 h) (cons 1 str))))

(defun dp:parse-areas (s)
  (setq s (vl-string-translate ",;\t" "   " s))
  (vl-remove-if-not 'numberp (read (strcat "(" s ")"))))

(defun dp:r3 (x) (rtos x 2 3))

(defun dp:shift (pts dx) (mapcar '(lambda (p) (list (+ (car p) dx) (cadr p))) pts))

;;; ------------------------------------------------------------------
;;;  MAIN COMMAND
;;; ------------------------------------------------------------------
(defun c:DIEPROFILE (/ *error* doc oldvars name sel shape lst s k bad
                       around afinal r rf n P G ap f scaled cent cc C
                       lay lay3 h stages i area pts e aact status rows ty
                       do3d ang xoff z prevreq req secs before solid)
  (setq doc (vla-get-ActiveDocument (vlax-get-acad-object)))
  (defun *error* (msg)
    (if (not (wcmatch (strcase msg) "*CANCEL*,*QUIT*,*EXIT*"))
      (princ (strcat "\nError: " msg)))
    (foreach v oldvars (setvar (car v) (cdr v)))
    (vla-EndUndoMark doc)
    (princ))
  (vla-StartUndoMark doc)
  (setq n *DP-NPTS*)

  ;; ---- 1. name ---------------------------------------------------
  (setq name (strcase (getstring "\nName for this die set (e.g. AIP) <DIE>: ")))
  (if (= name "") (setq name "DIE"))

  ;; ---- 2. final shape -------------------------------------------
  (while (null shape)
    (setq sel (entsel "\nSelect the FINAL shaped profile (closed polyline/spline): "))
    (cond
      ((null sel) (princ "\nNothing selected - try again (Esc to quit)."))
      ((not (member (cdr (assoc 0 (entget (car sel))))
                    '("LWPOLYLINE" "POLYLINE" "SPLINE" "ELLIPSE")))
       (princ "\nNot a polyline/spline. Use JOIN (or BOUNDARY) to make one closed polyline first."))
      ((not (vlax-curve-isClosed (car sel)))
       (princ "\nThat profile is not CLOSED. Close it (PEDIT > Close / JOIN) and try again."))
      (T (setq shape (car sel)))))

  ;; ---- 3. all areas, initial (round) -> final (shape) ------------
  (princ "\nEnter die areas in mm2 from INITIAL (round die) to FINAL (shape).")
  (princ "\nOne per line, or several separated by spaces. Empty Enter = finished.")
  (while (/= "" (setq s (getstring T (strcat "\n  Area of die " (itoa (1+ (length lst)))
                                             (if lst " <Enter = finish>: " ": ")))))
    (setq lst (append lst (mapcar 'float (dp:parse-areas s)))))
  (if (< (length lst) 2)
    (progn (princ "\nNeed at least 2 areas (initial round + final shape).") (exit)))
  (setq k 0)
  (while (< k (1- (length lst)))
    (if (>= (nth (1+ k) lst) (nth k lst)) (setq bad T))
    (setq k (1+ k)))
  (if bad
    (progn (princ "\nAreas must get SMALLER from initial to final. Please check and run again.") (exit)))
  (setq around (car lst)
        afinal (last lst)
        r      (sqrt (/ around pi)))
  (princ (strcat "\nRound die: area " (dp:r3 around) "  ->  diameter " (rtos (* 2 r) 2 3)))

  ;; ---- 4. shape points, final-area check -------------------------
  (setq P (dp:sample shape n))
  (if (< (dp:area P) 0.0) (setq P (reverse P)))
  (setq G  (dp:centroid P)
        ap (dp:area P))
  (princ (strcat "\nSelected shape area: " (dp:r3 (vlax-curve-getArea shape))
                 "   Final area entered: " (dp:r3 afinal)))
  (if (> (abs (- ap afinal)) (* afinal (/ *DP-TOL* 100.0)))
    (progn
      (princ "\n** The shape area does NOT match the final area.")
      (initget "Yes No")
      (if (/= (getkword "\nScale a copy of the shape to the final area? [Yes/No] <Yes>: ") "No")
        (setq f      (sqrt (/ afinal ap))
              P      (mapcar '(lambda (p) (list (+ (car G)  (* f (- (car p)  (car G))))
                                                (+ (cadr G) (* f (- (cadr p) (cadr G))))))
                             P)
              ap     (dp:area P)
              scaled T))))
  (if (<= around ap)
    (progn (princ "\nInitial (round) area must be bigger than the shape area.") (exit)))

  ;; ---- 5. centre and rounding ----------------------------------
  (setq cent (getpoint (strcat "\nPick die centre <Enter = shape centroid "
                               (rtos (car G) 2 3) "," (rtos (cadr G) 2 3) ">: ")))
  (setq cc (if cent (list (car cent) (cadr cent)) G))
  (initget 4)
  (setq rf (getreal "\nCorner rounding % (0 = no rounding) <6>: "))
  (if (null rf) (setq rf 6.0))
  (setq rf (/ rf 100.0))

  ;; ---- 6. build all stages (2D) ---------------------------------
  (princ "\nCalculating ...")
  (setq C    (dp:circpts P G cc r)
        lay  (strcat "DIE-" name)
        lay3 (strcat "DIE-" name "-3D")
        h    (/ r 10.0)
        i    1)
  ;; stages: list of (label area pts status) ; pts nil = round circle
  (foreach area lst
    (setq status nil pts nil)
    (cond
      ((= i 1)
       (dp:circle cc r lay 0.0)
       (setq status "round die (drawn)" aact (* pi r r)))
      ((= i (length lst))
       (setq pts P)
       (if scaled
         (setq e (dp:pline P lay 0.0) aact (vlax-curve-getArea e)
               status "final shape (scaled copy drawn)")
         (setq aact (vlax-curve-getArea shape) status "final shape (selected)")))
      ((<= area ap)
       (setq status "SKIPPED: not bigger than final shape" aact nil))
      (T
       (setq pts  (dp:profile P C area n rf)
             e    (dp:pline pts lay 0.0)
             aact (vlax-curve-getArea e)
             status "drawn")))
    (setq stages (cons (list (strcat name "-" (itoa i)) area pts status aact) stages)
          i (1+ i)))
  (setq stages (reverse stages))

  ;; ---- 7. labels + report ---------------------------------------
  (setq ty (+ (cadr cc) r))
  (dp:text (list (- (car cc) (* 1.5 h)) (+ (cadr cc) r (* 1.0 h))) (* 1.4 h) name lay)
  (princ (strcat "\n\n====== " name " ======"
                 "\n DIE NO        TARGET     DRAWN AREA    STATUS"))
  (foreach st stages
    (setq ty (- ty (* 1.8 h)))
    (dp:text (list (+ (car cc) r (* 2 h)) ty) h
             (strcat (car st) "  AREA " (dp:r3 (cadr st)) " SQ MM"
                     (if (wcmatch (cadddr st) "round*") (strcat "  ROUND DIA " (rtos (* 2 r) 2 3)) ""))
             lay)
    (princ (strcat "\n " (car st) "      " (dp:r3 (cadr st)) "     "
                   (if (nth 4 st) (rtos (nth 4 st) 2 4) "   -   ")
                   "     " (cadddr st))))

  (if (> (abs (- ap afinal)) (* afinal (/ *DP-TOL* 100.0)))
    (princ (strcat "\n** WARNING: final shape area " (dp:r3 ap) " is not the final area "
                   (dp:r3 afinal) " you entered. Check the shape size/scale!")))

  ;; ---- 8. optional 3D model ------------------------------------
  (initget "Yes No")
  (setq do3d (getkword "\nAlso create a separate 3D model of the process? [Yes/No] <Yes>: "))
  (if (/= do3d "No")
    (progn
      (initget 6)
      (setq ang (getreal "\nDie included angle in degrees (sets the length between dies) <16>: "))
      (if (null ang) (setq ang 16.0))
      (setq xoff (* 7.0 r) z 0.0 prevreq nil secs nil)
      (foreach st stages
        (if (wcmatch (cadddr st) "drawn,round*,final*")
          (progn
            (setq req (sqrt (/ (cadr st) pi)))          ; equivalent radius
            (if prevreq
              (setq z (+ z (/ (- prevreq req) (/ (sin (/ (* ang pi) 360.0))
                                                 (cos (/ (* ang pi) 360.0)))))))
            (setq prevreq req)
            (setq secs (cons (if (caddr st)
                               (dp:pline (dp:shift (caddr st) xoff) lay3 z)
                               (dp:circle (list (+ (car cc) xoff) (cadr cc)) r lay3 z))
                             secs)))))
      (setq secs (reverse secs))
      (dp:text (list (+ (car cc) xoff (- r)) (- (cadr cc) r (* 2.5 h))) (* 1.4 h)
               (strcat name " 3D  (" (rtos ang 2 1) " deg, length " (rtos z 2 2) ")") lay3)
      ;; loft through all sections, keep the section rings visible
      (setq oldvars (list (cons "CMDECHO" (getvar "CMDECHO"))
                          (cons "DELOBJ"  (getvar "DELOBJ"))
                          (cons "CLAYER"  (getvar "CLAYER"))))
      (setvar "CMDECHO" 0)
      (setvar "DELOBJ" 0)
      (setvar "CLAYER" lay3)
      (setq before (entlast))
      (command "_.LOFT")
      (foreach e secs (command e))
      (command "")
      (while (> (getvar "CMDACTIVE") 0) (command ""))
      (foreach v oldvars (setvar (car v) (cdr v)))
      (setq oldvars nil solid (entlast))
      (if (and solid (not (eq solid before)))
        (princ (strcat "\n3D solid created on layer " lay3 ", total length "
                       (rtos z 2 3) ".  View it: type  -VIEW  SWISO  then  VSCURRENT  S"))
        (princ (strcat "\n** LOFT did not succeed. The 3D sections are on layer " lay3
                       " - you can run LOFT on them manually.")))))

  (princ (strcat "\n2D profiles are on layer " lay ".  Type U to undo everything."))
  (vla-EndUndoMark doc)
  (princ))

(princ "\nDIEPROFILE loaded.  Type DIEPROFILE to run.")
(princ)
