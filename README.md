# DIEPROFILE – Wire-drawing die profile generator for AutoCAD
AutoLISP program that starts from a final shaped profile (Z / half-circle /
segment) and a list of die areas, and draws every intermediate die profile
from the round start die to the final shape, each with exactly its target
area. It can also build a 3D loft solid of the whole process.

## Use
1. In AutoCAD type `APPLOAD`, load `DIEPROFILE.lsp`.
2. Make the final shape one closed polyline (`JOIN` if needed).
3. Type `DIEPROFILE` and answer the prompts:
   - set name (any, e.g. AIP)
   - select the final shape
   - areas in mm² from INITIAL (round) to FINAL (shape); press Enter on an empty line to finish
   - die centre (Enter = shape centroid), corner rounding (Enter = 6 %)
   - 3D model? (Enter = Yes), die included angle (Enter = 16°)
4. Check the target vs. drawn areas printed on the command line.
   Type `U` to undo the whole run.

## How it works
Each point of the final shape moves radially toward the round die by a
fraction t. Corners are rounded progressively, and t is solved so that the
polyline area equals the target area (error < 0.000001 mm²). In the 3D
model, dies are spaced by Δ(equivalent radius) / tan(angle/2).
