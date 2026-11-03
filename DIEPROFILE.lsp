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
;;;                die sections a fixed distance apart along Z.
;;;
;;;  Method   : profiles are built from the final shape OUTWARD.  In every
;;;             direction from the die centre each profile moves toward the
;;;             round die by a fraction t, corners are rounded a little more
;;;             at each stage, and each profile is kept a small gap OUTSIDE
;;;             the previous one (no overlap).  t is solved so the drawn
;;;             polyline area equals the target (error < 0.000001 mm2).
;;;             AutoCAD's own area of every polyline is read back and printed.
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

;;; ------------------------------------------------------------------
;;;  POLAR method (used whenever possible) - guarantees NO OVERLAP
;;;  Every profile is described by its radius r(angle) from origin o,
;;;  at n equal angles.  Each profile is forced to stay outside the
;;;  previous (smaller) one by a small gap, so profiles never touch or
;;;  cross, and the area is then solved exactly with that rule active.
;;; ------------------------------------------------------------------

;; Radius of a closed CCW polygon at n equal angles from origin o.
;; Returns (a0 r0 r1 ... r[n-1]) or nil if some direction from o
;; crosses the outline more than once ("not star-shaped" from o).
(defun dp:polar-shape (pts o n / lst ok sum s a0 d segs p q k tau w u den r res)
  (setq lst (append pts (list (car pts))) ok T sum 0.0
        s (angle o (car pts)) a0 s)
  (while (and ok (cdr lst))
    (setq p (car lst) q (cadr lst)
          d (dp:normang (- (angle o q) (angle o p))))
    (if (< d -1e-9) (setq ok nil))
    (setq segs (cons (list s (+ s d) p q) segs)
          s (+ s d) sum (+ sum d) lst (cdr lst)))
  (if (and ok (< (abs (- sum (* 2.0 pi))) 1e-6))
    (progn
      (setq segs (reverse segs) k 0)
      (repeat n
        (setq tau (+ a0 (* k (/ (* 2.0 pi) n))))
        (while (and (cdr segs) (> tau (cadr (car segs)))) (setq segs (cdr segs)))
        (setq p   (caddr (car segs)) q (cadddr (car segs))
              w   (list (- (car q) (car p)) (- (cadr q) (cadr p)))
              u   (list (cos tau) (sin tau))
              den (- (* (car u) (cadr w)) (* (cadr u) (car w)))
              r   (if (< (abs den) 1e-12)
                    (distance o p)
                    (/ (- (* (- (car p) (car o)) (cadr w)) (* (- (cadr p) (cadr o)) (car w)))
                       den))
              res (cons r res) k (1+ k)))
      (cons a0 (reverse res)))))

;; Same, for ANY shape (slots, hooks): radius of the OUTERMOST crossing in
;; each direction.  The region inside it always contains the whole shape.
(defun dp:polar-env (pts o n / lst segs p q w po k tau u den s r best res)
  (setq lst (append pts (list (car pts))))
  (while (cdr lst)
    (setq p (car lst) q (cadr lst)
          segs (cons (list (list (- (car q) (car p)) (- (cadr q) (cadr p)))
                           (list (- (car p) (car o)) (- (cadr p) (cadr o))))
                     segs)
          lst (cdr lst)))
  (setq k 0)
  (repeat n
    (setq tau (* k (/ (* 2.0 pi) n)) u (list (cos tau) (sin tau)) best 0.0)
    (foreach sg segs
      (setq w (car sg) po (cadr sg)
            den (- (* (car u) (cadr w)) (* (cadr u) (car w))))
      (if (> (abs den) 1e-12)
        (progn
          (setq s (/ (- (* (car po) (cadr u)) (* (cadr po) (car u))) den)
                r (/ (- (* (car po) (cadr w)) (* (cadr po) (car w))) den))
          (if (and (>= s -1e-9) (<= s (+ 1.0 1e-9)) (> r best)) (setq best r)))))
    (setq res (cons best res) k (1+ k)))
  (cons 0.0 (reverse res)))

;; distance from o to the round-die circle (centre cc, radius r) per angle
(defun dp:Rpolar (o cc r n a0 / v vv k tau b res)
  (setq v  (list (- (car cc) (car o)) (- (cadr cc) (cadr o)))
        vv (+ (* (car v) (car v)) (* (cadr v) (cadr v))) k 0)
  (repeat n
    (setq tau (+ a0 (* k (/ (* 2.0 pi) n)))
          b   (+ (* (cos tau) (car v)) (* (sin tau) (cadr v)))
          res (cons (+ b (sqrt (max 0.0 (- (+ (* b b) (* r r)) vv)))) res)
          k   (1+ k)))
  (reverse res))

;; exact area of the polygon through the n polar points
(defun dp:parea (rs n / l s)
  (setq l (append rs (list (car rs))) s 0.0)
  (while (cdr l) (setq s (+ s (* (car l) (cadr l))) l (cdr l)))
  (* 0.5 (sin (/ (* 2.0 pi) n)) s))

;; closed moving average of a list of numbers, 3 passes
(defun dp:box1 (l m n / w ext lead tail sm res)
  (setq w (+ m m 1) ext (append (dp:lastn l m) l (dp:firstn l m))
        lead ext tail ext sm 0.0)
  (repeat w (setq sm (+ sm (car lead)) lead (cdr lead)))
  (repeat n
    (setq res (cons (/ sm w) res))
    (if lead (setq sm (+ sm (- (car lead) (car tail))) lead (cdr lead)))
    (setq tail (cdr tail)))
  (reverse res))

(defun dp:smooth1 (l m n)
  (if (>= m 1) (repeat 3 (setq l (dp:box1 l m n))))
  l)

;; profile radii at fraction tt:  between floor fl and the circle Rl
(defun dp:pp (S1 S2 Rl fl tt)
  (mapcar '(lambda (a b c f) (min c (max (+ a (* tt (- b a))) f))) S1 S2 Rl fl))

(defun dp:psolve (S1 S2 Rl fl target n / lo hi mid)
  (setq lo 0.0 hi 1.0)
  (repeat 50
    (setq mid (/ (+ lo hi) 2.0))
    (if (< (dp:parea (dp:pp S1 S2 Rl fl mid) n) target) (setq lo mid) (setq hi mid)))
  (/ (+ lo hi) 2.0))

;; one stage: shape radii rs, circle radii Rl, previous (inner) profile
;; prev with area prevA  ->  radii of new profile with area = target
(defun dp:pstage (rs Rl prev prevA target n rf / gap fl k t0 m S1 S2)
  (setq gap (* 0.25 (- (sqrt (/ target pi)) (sqrt (/ prevA pi)))) k 0)
  (setq fl (mapcar '(lambda (a c) (min c (+ a gap))) prev Rl))
  (while (and (< k 30) (>= (dp:parea fl n) target))     ; gap too big -> reduce
    (setq gap (/ gap 2.0) k (1+ k)
          fl  (mapcar '(lambda (a c) (min c (+ a gap))) prev Rl)))
  (setq t0 (dp:psolve rs Rl Rl fl target n)
        m  (fix (+ 0.5 (* t0 n rf))))
  (if (> m (/ n 4)) (setq m (/ n 4)))
  (setq S1 (dp:smooth1 rs m n) S2 (dp:smooth1 Rl m n))
  (dp:pp S1 S2 Rl fl (dp:psolve S1 S2 Rl fl target n)))

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
                       pol tmp o a0 rs Rl angs prev prevA mids rr envA
                       do3d dz xoff z secs before solid)
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
  (setq lay  (strcat "DIE-" name)
        lay3 (strcat "DIE-" name "-3D")
        h    (/ r 10.0))
  ;; polar description of the shape, seen from die centre (or centroid)
  (foreach cand (list cc G)
    (if (and (null pol) (setq tmp (dp:polar-shape P cand n)))
      (setq pol tmp o cand)))
  (if (null pol)
    (progn
      (setq pol (dp:polar-env P cc n) o cc)
      (princ (strcat "\n** NOTE: the shape has a slot/hook (some directions from the centre"
                     "\n   cross its edge twice). The first profile will bridge over it."))))
  (setq a0 (car pol) rs (cdr pol) Rl (dp:Rpolar o cc r n a0) k 0 angs nil)
  (repeat n (setq angs (cons (+ a0 (* k (/ (* 2.0 pi) n))) angs) k (1+ k)))
  (setq angs (reverse angs) envA (dp:parea rs n))
  (if (vl-some '(lambda (a c) (> a c)) rs Rl)
    (princ "\n** WARNING: part of the shape is OUTSIDE the round die circle!"))
  ;; intermediate profiles, built from the final shape OUTWARD so each
  ;; one is kept outside the previous one
  (setq prev rs prevA envA mids nil i (1- (length lst)))
  (while (>= i 2)
    (setq area (nth (1- i) lst))
    (if (> area prevA)
      (setq rr    (dp:pstage rs Rl prev prevA area n rf)
            prev  rr
            prevA area
            mids  (cons (cons i (mapcar '(lambda (x a) (polar o a x)) rr angs)) mids)))
    (setq i (1- i)))
  (setq i 1)
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
      ((null (setq pts (cdr (assoc i mids))))
       (setq status "SKIPPED: not bigger than final shape" aact nil))
      (T
       (setq e    (dp:pline pts lay 0.0)
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
      (setq dz (getreal "\nDistance between die sections in 3D <5>: "))
      (if (null dz) (setq dz 5.0))
      (setq xoff (* 7.0 r) z nil secs nil)
      (foreach st stages
        (if (wcmatch (cadddr st) "drawn,round*,final*")
          (progn
            (setq z (if z (+ z dz) 0.0))
            (setq secs (cons (if (caddr st)
                               (dp:pline (dp:shift (caddr st) xoff) lay3 z)
                               (dp:circle (list (+ (car cc) xoff) (cadr cc)) r lay3 z))
                             secs)))))
      (setq secs (reverse secs))
      (dp:text (list (+ (car cc) xoff (- r)) (- (cadr cc) r (* 2.5 h))) (* 1.4 h)
               (strcat name " 3D  (section spacing " (rtos dz 2 2) ")") lay3)
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
