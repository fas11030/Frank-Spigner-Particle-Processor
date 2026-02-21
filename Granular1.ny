;nyquist plug-in
;version 4
;type process
;name "Frank Spigner Particle Processor)"
;action "Processing..."
;info "Version notes: Output Mode maps 1:1 to design diagram blocks.\nUses timed-seq (score) to avoid sim/simrep stack overflows.\nOutput duration ALWAYS equals the selection (process-safe)."

;control mode "Stage" choice "Final,Selected,DensityClicks,TimeMapClicks,Extract,Envelope,Amplitude,TimeTransform,FreqMult,FreqLattice" 0

;control density  "Density Multiplier (x)"              real "" 1.0   0.01 40.0
;control durmult  "Duration Multiplier (source map x)"  real "" 1.0   0.10 10.0
;control maxgr    "Max Particles (hard cap)"            int  "" 2000  1 200000

;control gmin     "Particle MIN Length (samples)"       int  "" 4000  4 500000
;control gmax     "Particle MAX Length (samples)"       int  "" 12000 4 500000
;control gdist    "Particle Length Distribution"        choice "TruncGaussian,Uniform,Triangular,Cauchy,WignerSemicircle" 1

;control tjmax    "Time Transform Range (± samples)"    int  "" 0     0 200000
;control tjdist   "Time Transform Distribution"         choice "TruncGaussian,Uniform,Triangular,Cauchy,WignerSemicircle" 1

;control amin     "Attack MIN (samples)"                int "" 50    1 500000
;control amax     "Attack MAX (samples)"                int "" 400   1 500000
;control dlmin    "Delay  MIN (samples)"                int "" 1     1 500000
;control dlmax    "Delay  MAX (samples)"                int "" 200   1 500000
;control dcmin    "Decay  MIN (samples)"                int "" 20    1 500000
;control dcmax    "Decay  MAX (samples)"                int "" 400   1 500000
;control smin     "Sustain MIN (samples)"               int "" 1     1 500000
;control smax     "Sustain MAX (samples)"               int "" 2000  1 500000
;control rmin     "Release MIN (samples)"               int "" 200   1 500000
;control rmax     "Release MAX (samples)"               int "" 2000  1 500000
;control slevel   "Sustain Level (0-1)"                 real "" 0.8  0.0 1.0
;control edist    "Envelope Stage Distribution"         choice "TruncGaussian,Uniform,Triangular,Cauchy,WignerSemicircle" 1

;control gdbmin   "Particle Gain MIN (dB)"              real "" -6.0  -60.0 24.0
;control gdbmax   "Particle Gain MAX (dB)"              real ""  0.0  -60.0 24.0
;control gdbdist  "Particle Gain Distribution"          choice "TruncGaussian,Uniform,Triangular,Cauchy,WignerSemicircle" 1

;control fminmul  "Freq Mult MIN (x)"                   real "" 1.0   0.25 8.0
;control fmaxmul  "Freq Mult MAX (x)"                   real "" 1.0   0.25 8.0
;control fdist    "Freq Mult Distribution"              choice "TruncGaussian,Uniform,Triangular,Cauchy,WignerSemicircle" 1

;control fund     "Lattice Fundamental (Hz)"            real "Hz" 440.0 10.0 20000.0
;control order    "Max Harm/Subharm Order"              int  "" 8     1 64
;control nodefmin "Allowed Node MIN (Hz)"               real "Hz" 20.0   1.0 20000.0
;control nodefmax "Allowed Node MAX (Hz)"               real "Hz" 20000.0 1.0 20000.0
;control spreadhz "Node Spread (±Hz)"                   real "Hz" 0.0    0.0 5000.0
;control ndist    "Node Choice Distribution"            choice "TruncGaussian,Uniform,Triangular,Cauchy,WignerSemicircle" 1
;control sdist    "Node Spread Distribution"            choice "TruncGaussian,Uniform,Triangular,Cauchy,WignerSemicircle" 1

;control outgain  "Output Gain (dB)"                    real "" 0.0  -36.0 24.0


;; ---------------- numeric helpers ----------------
(setq *pi* 3.141592653589793)

(defun clamp (x lo hi) (max lo (min hi x)))
(defun mag (x) (if (< x 0.0) (- x) x))
(defun u01 () (rrandom)) ; 0..1

(defun gauss0 ()
  (let* ((u1 (max 1.0e-12 (u01)))
         (u2 (u01)))
    (* (sqrt (* -2.0 (log u1)))
       (cos (* 2.0 *pi* u2)))))

(defun disk-x ()
  (let* ((x (- (* 2.0 (u01)) 1.0))
         (y (- (* 2.0 (u01)) 1.0)))
    (if (<= (+ (* x x) (* y y)) 1.0) x (disk-x))))

(defun dist-sym (dtype)
  ;; returns in [-1,1]
  (cond
    ((= dtype 0) ; TruncGaussian
     (let ((x 2.0))
       (while (> (mag x) 1.0) (setq x (gauss0)))
       x))
    ((= dtype 1) ; Uniform
     (- (* 2.0 (u01)) 1.0))
    ((= dtype 2) ; Triangular
     (- (+ (u01) (u01)) 1.0))
    ((= dtype 3) ; Cauchy (truncated)
     (let ((x 2.0))
       (while (> (mag x) 1.0)
         (setq x (tan (* *pi* (- (u01) 0.5)))))
       x))
    (t           ; Wigner semicircle
     (disk-x))))

(defun sample-range (dtype lo hi)
  (let* ((a (min lo hi)) (b (max lo hi)))
    (+ a (* 0.5 (+ 1.0 (dist-sym dtype)) (- b a)))))

(defun sample-int (dtype lo hi)
  (round (sample-range dtype lo hi)))

(defun wrap0 (x period)
  ;; wrap into [0,period) without FLOOR
  (if (<= period 0.0)
      0.0
      (let ((y x))
        (while (< y 0.0) (setq y (+ y period)))
        (while (>= y period) (setq y (- y period)))
        y)))

(defun get-sel-samps (sig sr)
  ;; Audacity provides LEN = selection samples in Effects context.
  (cond
    ((and (boundp 'len) (numberp len) (> len 0)) (truncate len))
    (t
     (let ((st (snd-stop-time sig)))
       (if (and (numberp st) (> st 0.0) (< st 1.0e20))
           (truncate (* st sr))
           (error "Cannot determine selection length (LEN not set and stop-time not finite)."))))))


;; ---------------- early blocks: click ----------------
(defun click-sound (src sr srcdur)
  ;; 1ms excerpt from mid selection
  (let* ((cdur (max (/ 4.0 sr) 0.001))
         (mid (max 0.0 (- (* 0.5 srcdur) (* 0.5 cdur))))
         (maxp (max 0.0 (- srcdur cdur)))
         (p (clamp mid 0.0 maxp)))
    (cue (extract-abs p (+ p cdur) src))))


;; ---------------- time mapping blocks ----------------
(defun mapped-srcpos (outpos gdur sr srcdur use-jitter)
  ;; Block 2: outpos -> source position via durmult
  ;; Block 6: add jitter (±tjmax samples)
  (let* ((base (/ outpos (max 0.001 durmult)))
         (jit  (if use-jitter
                   (/ (round (* tjmax (dist-sym tjdist))) sr)
                   0.0))
         (x (+ base jit))
         (maxsrc (max 0.0 (- srcdur gdur))))
    (if (<= maxsrc 0.0) 0.0 (wrap0 x maxsrc))))


;; ---------------- envelope block (true A + Delay + Decay + Sustain + Release) ----------------
(defun pick-adsdrs (glen)
  ;; pick A, DL, DC, S, R in samples; then scale down until sum <= glen
  (let* ((a  (max 1 (sample-int edist amin  amax)))
         (dl (max 1 (sample-int edist dlmin dlmax)))
         (dc (max 1 (sample-int edist dcmin dcmax)))
         (s  (max 1 (sample-int edist smin  smax)))
         (r  (max 1 (sample-int edist rmin  rmax)))
         (sum (+ a dl dc s r)))
    (if (> sum glen)
        (let* ((fac (/ (max 5 glen) (max 1.0 sum))))
          (setq a  (max 1 (round (* a  fac))))
          (setq dl (max 1 (round (* dl fac))))
          (setq dc (max 1 (round (* dc fac))))
          (setq s  (max 1 (round (* s  fac))))
          (setq r  (max 1 (round (* r  fac))))
          ;; still too long? steal from sustain first
          (while (> (+ a dl dc s r) glen)
            (if (> s 1) (setq s (- s 1))
                (if (> dl 1) (setq dl (- dl 1))
                    (if (> dc 1) (setq dc (- dc 1))
                        (if (> r 1) (setq r (- r 1))
                            (if (> a 1) (setq a (- a 1))))))))))
    ;; if too short due to rounding, pad sustain
    (while (< (+ a dl dc s r) glen) (setq s (+ s 1)))
    (list a dl dc s r)))

(defun adsdrs-env (sr glen)
  ;; Generate envelope at AUDIO SR to avoid FORCE-SRATES noise
  (let* ((lst (pick-adsdrs glen))
         (aS  (nth 0 lst))
         (dlS (nth 1 lst))
         (dcS (nth 2 lst))
         (sS  (nth 3 lst))
         (rS  (nth 4 lst))
         (aT  (/ aS sr))
         (dlT (/ dlS sr))
         (dcT (/ dcS sr))
         (sT  (/ sS sr))
         (rT  (/ rS sr))
         (t0  0.0)
         (tA  (+ t0 aT))
         (tDL (+ tA dlT))
         (tDC (+ tDL dcT))
         (tS  (+ tDC sT))
         (tE  (+ tS rT))
         (sv  (clamp slevel 0.0 1.0)))
    (control-srate-abs sr
      (let ((*warp* 1.0))
        (pwl t0 0.0
             tA 1.0
             tDL 1.0
             tDC sv
             tS  sv
             tE  0.0)))))


;; ---------------- amplitude block ----------------
(defun pick-gain-db () (sample-range gdbdist gdbmin gdbmax))


;; ---------------- frequency blocks ----------------
(defun pick-fmul () (max 0.001 (sample-range fdist fminmul fmaxmul)))

(defun pick-node-hz ()
  (let* ((ord (max 1 order))
         (k (max 1 (min ord (round (sample-range ndist 1.0 ord)))))
         (harm? (< (u01) 0.5))
         (base (if harm? (* fund k) (/ fund k)))
         (dHz (if (> spreadhz 0.0) (sample-range sdist (- spreadhz) spreadhz) 0.0))
         (hz (+ base dHz))
         (lo (min nodefmin nodefmax))
         (hi (max nodefmin nodefmax)))
    (clamp hz lo hi)))

(defun pitch-grain (g sr ratio gdur)
  ;; fixed-ratio pitch via scale-srate + resample; keep duration by pad/trim
  (let* ((ratio (max 0.001 ratio)))
    (if (= ratio 1.0)
        g
        (let* ((pg (resample (scale-srate g ratio) sr))
               (pg (cue pg))
               (pd (snd-stop-time pg)))
          (cond
            ((< pd gdur) (cue (seq pg (s-rest (- gdur pd)))))
            ((> pd gdur) (cue (extract-abs 0.0 gdur pg)))
            (t pg))))))


;; ---------------- particle builder (stages 3..8) ----------------
(defun build-particle (src sr srcdur outpos outdur glen stage)
  ;; stage:
  ;; 4 Extract
  ;; 5 + Envelope
  ;; 6 + Gain
  ;; 7 + Time Transform (jitter)
  ;; 8 + Freq Mult
  ;; 9 + Freq Lattice
  (let* ((glen (max 4 glen))
         (gdur (/ glen sr))
         ;; keep grain within remaining output
         (rem (max 0.0 (- outdur outpos)))
         (gdur (min gdur rem))
         (gdur (max (/ 4.0 sr) gdur))
         (glen (max 4 (round (* gdur sr))))

         (use-jitter (>= stage 7))
         (srcpos (mapped-srcpos outpos gdur sr srcdur use-jitter))

         (g (cue (extract-abs srcpos (+ srcpos gdur) src)))

         ;; freq blocks (approx by pitch ratio)
         (fmul (pick-fmul))
         (nodehz (pick-node-hz))
         (ratio (cond
                  ((>= stage 9) (* fmul (/ nodehz (max 1.0e-9 fund))))
                  ((>= stage 8) fmul)
                  (t 1.0)))
         (g (if (>= stage 8) (pitch-grain g sr ratio gdur) g))

         ;; env + gain
         (e (adsdrs-env sr glen))
         (g (if (>= stage 5) (mult g e) g))

         (gdB (pick-gain-db))
         (g (if (>= stage 6) (scale-db gdB g) g)))
    g))


;; ---------------- timed-seq mixer (FIX: behaviors must be expressions like (cue <sound>)) ----------------
(defun make-event (t snd)
  ;; TIMED-SEQ requires each event = (time stretch behavior).
  ;; The behavior must respond to *warp*, so we use (cue <sound>), not a bare sound object.
  (list t 1.0 (list 'cue snd)))

(defun render-score (bed events outdur)
  (let* ((score (reverse events))
         (mix (timed-seq score)))
    (cue (extract-abs 0.0 outdur (sum bed mix)))))


;; ---------------- per-channel processor ----------------
(defun process-one (raw)
  (let* ((sig   (cue raw))
         (sr    (snd-srate sig))
         (ns    (max 4 (get-sel-samps sig sr)))
         (srcdur (/ ns sr))
         (src   (cue (extract-abs 0.0 srcdur sig)))

         ;; process-safe output
         (outdur srcdur)
         (bed (mult 0.0 src))

         (gminS (clamp gmin 4 ns))
         (gmaxS (clamp gmax gminS ns))

         (avgS  (* 0.5 (+ gminS gmaxS)))
         (avgDur (/ avgS sr))
         (hop0  (max (/ 4.0 sr) (/ avgDur (max 0.001 density))))
         (wantN (max 1 (truncate (/ outdur hop0))))
         (ngr   (min maxgr wantN))
         (hop   (/ outdur ngr))

         (clk (click-sound src sr srcdur))
         (events nil))

    ;; Stage 1: selected audio (no score)
    (if (= mode 1)
        (cue (scale-db outgain src))

        (progn
          (dotimes (i ngr)
            (let* ((outpos (* i hop))
                   (remS (max 4 (truncate (* sr (max 0.0 (- outdur outpos))))))
                   (glen0 (clamp (sample-int gdist gminS gmaxS) 4 ns))
                   (glen  (min glen0 remS))
                   (snd nil))

              (cond
                ;; Stage 2: density clicks (no mapping)
                ((= mode 2)
                 (setq snd clk))

                ;; Stage 3: time-map clicks (mapped source position, no jitter)
                ((= mode 3)
                 (let* ((cdur (snd-stop-time clk))
                        (srcp (mapped-srcpos outpos cdur sr srcdur nil))
                        (sndc (cue (extract-abs srcp (+ srcp cdur) src))))
                   (setq snd sndc)))

                ;; Stage 4..9: cumulative grain stages
                ((= mode 4) (setq snd (build-particle src sr srcdur outpos outdur glen 4)))
                ((= mode 5) (setq snd (build-particle src sr srcdur outpos outdur glen 5)))
                ((= mode 6) (setq snd (build-particle src sr srcdur outpos outdur glen 6)))
                ((= mode 7) (setq snd (build-particle src sr srcdur outpos outdur glen 7)))
                ((= mode 8) (setq snd (build-particle src sr srcdur outpos outdur glen 8)))
                ((= mode 9) (setq snd (build-particle src sr srcdur outpos outdur glen 9)))

                ;; Stage 0: Final mix (all blocks)
                (t (setq snd (build-particle src sr srcdur outpos outdur glen 9))))

              (setq events (cons (make-event outpos snd) events))))

          (scale-db outgain (render-score bed events outdur))))))

;; ---------------- input handling (s / *track*) ----------------
(setq in
  (cond
    ((and (boundp 's) (or (soundp s) (arrayp s))) s)
    ((and (boundp '*track*) (or (soundp *track*) (arrayp *track*))) *track*)
    (t (error "No input sound found"))))

(if (arrayp in)
    (let* ((n (length in))
           (outarr (make-array n)))
      (dotimes (k n outarr)
        (setf (aref outarr k) (process-one (aref in k)))))
    (process-one in))