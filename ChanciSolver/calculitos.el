;; This buffer is for notes you don't want to save, and for Lisp evaluation.
;; If you want to create a file, visit that file with C-x C-f,
;; then enter the text in that file's own buffer.

(+ 2 2)

(defun bo(pressure)
  (/ 1 (+ 1 (* 0.000025 (- pressure 3000)))))

(defun convpress (pressure)
  (* 6894.76 pressure))

(mapcar #'bo (list 2900 2500 2000 1500 1000))
(mapcar #'convpress (list 2900 2500 2000 1500 1000))


(1.0025062656641603 1.0126582278481011 1.0256410256410258 1.0389610389610389 1.0526315789473684)
(19994804.0 17236900.0 13789520.0 10342140.0 6894760.0)

(convpress 3000)
