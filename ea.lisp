;;;; Copyright (c) 2012, Christopher Mark Gore,
;;;; All rights reserved.
;;;;
;;;; 8729 Lower Marine Road, Saint Jacob, Illinois 62281 USA.
;;;; Web: http://www.cgore.com
;;;; Email: cgore@cgore.com
;;;;
;;;; Redistribution and use in source and binary forms, with or without
;;;; modification, are permitted provided that the following conditions are met:
;;;;
;;;;     * Redistributions of source code must retain the above copyright
;;;;       notice, this list of conditions and the following disclaimer.
;;;;
;;;;     * Redistributions in binary form must reproduce the above copyright
;;;;       notice, this list of conditions and the following disclaimer in the
;;;;       documentation and/or other materials provided with the distribution.
;;;;
;;;;     * Neither the name of Christopher Mark Gore nor the names of other
;;;;       contributors may be used to endorse or promote products derived from
;;;;       this software without specific prior written permission.
;;;;
;;;; THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
;;;; AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
;;;; IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
;;;; ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
;;;; LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
;;;; CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
;;;; SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
;;;; INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
;;;; CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
;;;; ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
;;;; POSSIBILITY OF SUCH DAMAGE.

(defvar *noisy* t)
(defun noisy (&rest rest)
  (if *noisy*
    (progn (apply #'format t rest)
           (sleep 0.25))))

(defun the-last (list)
  (assert (listp list))
  (car (last list)))

(defun random-in-range (lower upper)
  "(random-in-range lower upper) --> result.
  Produces a random number between the lower and upper bounds specified."
  (assert (numberp lower))
  (assert (numberp upper))
  (let ((a (min lower upper))
        (z (max lower upper)))
    (+ a (random (- z a)))))

(defun random-element (list)
  (assert (listp list))
  (nth (random (1- (length list))) list))

;;; Cf. http://www.necessaryandsufficient.net/2010/12/
;;;     classical-test-functions-for-genetic-algorithms/
;;;
;;; Cf. http://en.wikipedia.org/wiki/Rastrigin_function
(defvar *rastrigin-lower* -5.12)
(defvar *rastrigin-upper* 5.12)
(defun rastrigin (a x)
  (assert (numberp a)) ; The coefficient A should be a number.
  (assert (listp x)) ; The argument x should be a list.
  (mapcar (lambda (xi) ; Every x_i should be a number, correctly bounded.
	    (assert (numberp xi))
	    (assert (<= *rastrigin-lower* xi *rastrigin-upper*)))
          x)
  (let ((n (length x)))
    (+ (* a n)
       (reduce #'+ (mapcar (lambda (xi)
			     (- (expt xi 2.0)
				(* a (cos (* 2 pi xi)))))
			   x)))))

;; Some simple examples/unit tests of the Rastrigin function.
(assert (= 0.0 (rastrigin 10 '(0 0)))) ; This is the global minimum.
(assert (= 0.0 (rastrigin 10 '(0.0 0.0)))) ; This is the global minimum.
(assert (= 0.0 (rastrigin 10.0 '(0.0 0.0)))) ; This is the global minimum.
(assert (= 0.0 (rastrigin 10 '(0 0 0 0 0 0)))) ; This is the global minimum.

(defun bounded (min value max)
  "(bounded min value max) --> result.
  This is a simple function to bound a value between a minumum value and a
  maximum value.  If the value is within the bounds, then you get it back.
  If it is above the maximum, then you get the maximum instead.
  If it is below the minimum, then you get the minimum instead."
  (let ((a (min min max))
	(z (max min max)))
    (min (max a value) z)))

;; Some simple examples/unit tests of the bounded function.
(assert (= 55 (bounded 10 55 100))) ; Within the bounds.
(assert (= 10 (bounded 10 5 100))) ; Below the minimum.
(assert (= 100 (bounded 10 500 100))) ; Above the maximum.
(assert (= 10 (bounded 100 5 10))) ; Below the minimum, min/max swapped.
(assert (= 100 (bounded 100 500 10))) ; Above the maximum, min/max swapped.

(defclass gene () ; A single floating-point gene.
  ((value :accessor value :initarg :value
	  :type float :initform nil)
   (lower-bound :accessor lower-bound :initarg :lower-bound
		:type float :initform nil)
   (upper-bound :accessor upper-bound :initarg :upper-bound
		:type float :initform nil)
   (individual :accessor individual :initarg :individual
	       :type individual :initform nil)))

(defclass individual () ; A simple evolvable individual.
  ((genotype :accessor genotype :initarg :genotype
	     :type (vector gene) :initform nil)
   (generation :accessor generation :initarg :generation
               :type (integer 0 *) :initform nil)
   (fitness :accessor fitness :initarg :fitness
	    :type float :initform nil)
   (ea :accessor ea :initarg :ea
       :type ea :initform nil)))

(defclass ea () ; The base class for all of our EAs.
  ((generation :accessor generation :initarg :generation
               :type (integer 0 *) :initform 0)
   (population :accessor population :initarg :population
               :type list :initform nil)
   (parents :accessor parents :initarg :parents
            :type list :initform nil)
   (children :accessor children :initarg :children
             :type list :initform nil)
   (environ-lower :accessor environ-lower :initarg :environ-lower
                  :type list :initform nil)
   (environ-upper :accessor environ-upper :initarg :environ-upper
		  :type list :initform nil)
   (min-pop-size :accessor min-pop-size :initarg :min-pop-size
		 :type (integer 1 *) :initform 50)
   (max-pop-size :accessor max-pop-size :initarg :max-pop-size
		 :type (integer 1 *) :initform 50)
   (mutation-prob :accessor mutation-prob :initarg :mutation-prob
		  :type (float 0.0 1.0) :initform 0.1)
   (mutation-factor :accessor mutation-factor :initarg :mutation-factor
		    :type (float 0.0 1.0) :initform 0.05)
   (crossover-prob :accessor crossover-prob :initarg :crossover-prob
		   :type (float 0.0 1.0) :initform 0.2)))

(defclass easy2d-ea (ea) ; An easy problem EA.
  ((environ-lower :initform '(-10.0 -10.0))
   (environ-upper :initform '(10.0 10.0))))


(defclass rastrigin2d-ea (ea) ; An EA for the 2D Rastrigin function.
  ((environ-lower :initform (list *rastrigin-lower* *rastrigin-lower*))
   (environ-upper :initform (list *rastrigin-upper* *rastrigin-upper*))
   (coefficient-a :accessor coefficient-a :initarg :coefficient-a
                  :type float :initform 10.0)))

(defgeneric duplicate (thing))
(defgeneric mutate? (individual))
(defgeneric crossover? (individual))
(defgeneric mutate (thing))
(defgeneric mutate! (thing))
(defgeneric crossover (father mother))
(defgeneric make-random-individual (ea))
(defgeneric add-random-individual (ea))
(defgeneric fitness-function (individual ea))
(defgeneric initialize-random-population (ea))
(defgeneric evaluate-fitness (thing))
(defgeneric terminate-evolution? (ea))
(defgeneric select-parents (ea))
(defgeneric recombine-parents (ea))
(defgeneric mutate-children (ea))
(defgeneric integrate-children (ea))
(defgeneric select-survivors (ea))
(defgeneric evolve (ea))

(defmethod print-object ((gene gene) stream)
  (format stream "~A" (value gene)))

(defmethod print-object ((individual individual) stream)
  (format stream "<Individual: Generation ~A // Fitness ~A // "
          (generation individual) (fitness individual))
  (if (genotype individual)
    (loop for i from 0 to (1- (length (genotype individual))) do
          (format stream "~A " (aref (genotype individual) i)))
    "no genes")
  (format stream ">"))

(defmethod duplicate ((list list))
  (mapcar #'duplicate list))

(defmethod duplicate ((vector vector))
  (coerce (duplicate (coerce vector 'list)) 'vector))

(defmethod duplicate ((gene gene))
  (with-slots (value lower-bound upper-bound individual) gene
    (make-instance 'gene
		   :value value
                   :lower-bound lower-bound :upper-bound upper-bound
		   :individual individual)))

(defmethod duplicate ((individual individual))
  (with-slots (genotype fitness ea generation) individual
    (make-instance 'individual
                   :genotype (duplicate genotype) :fitness fitness :ea ea
                   :generation generation)))

(defmethod mutate? ((individual individual))
  (< (random 1.0) (mutation-prob (ea individual))))

(defmethod crossover? ((individual individual))
  (< (random 1.0) (crossover-prob (ea individual))))

(defmethod mutate! ((gene gene))
  (with-slots (value lower-bound upper-bound individual) gene
    (let* ((range (abs (- upper-bound lower-bound)))
	   (abs-delta (* range (random (/ (mutation-factor (ea individual))))))
	   (sign (if (= 0 (random 1)) -1 +1))
	   (delta (* sign abs-delta))
	   (new-value (+ value delta)))
      (setf value (bounded lower-bound new-value upper-bound))))
  gene)

(defmethod mutate ((gene gene))
  (mutate! (duplicate gene)))

(defmethod mutate! ((individual individual))
  (with-slots (genotype) individual
    (when (mutate? individual)
      (loop for i from 0 to (1- (length genotype)) do
	   (mutate! (aref genotype i)))))
  individual)

(defmethod mutate ((individual individual))
  (mutate! (duplicate individual)))

(defmethod crossover ((father gene) (mother gene))
  (let ((result (duplicate father)))
    (when (= 0 (random 1))
      (setf (value result) (value mother)))
    result))

(defmethod crossover ((father individual) (mother individual))
  (let ((result (duplicate father)))
    (setf (generation result) (generation (ea result)))
    (loop for i from 0 to (1- (length (genotype father))) do
	 (setf (aref (genotype result) i)
	       (crossover (aref (genotype father) i)
                          (aref (genotype mother) i))))
    result))

(defmethod make-random-individual ((ea ea))
  (let ((result (make-instance 'individual :ea ea :generation (generation ea))))
    (setf (genotype result)
          (coerce (loop for i from 0 to (1- (length (environ-lower ea))) collect
                        (make-instance 'gene 
                                       :value (random-in-range
                                                (nth i (environ-lower ea))
                                                (nth i (environ-upper ea)))
                                       :individual result
                                       :lower-bound (nth i (environ-lower ea))
                                       :upper-bound (nth i (environ-upper ea))))
                  'vector))
    result))

(defmethod add-random-individual ((ea ea))
  (setf (population ea)
	(cons (make-random-individual ea) (population ea)))
  (noisy "Adding a new random individual: ~A.~%" (first (population ea))))

(defmethod initialize-random-population ((ea ea))
  (setf (population ea) nil)
  (loop while (< (length (population ea)) (min-pop-size ea)) do
        (add-random-individual ea)))

(defmethod evaluate-fitness ((individual individual))
  (when (null (fitness individual))
    (setf (fitness individual)
          (fitness-function individual (ea individual)))))

(defmethod evaluate-fitness ((ea ea))
  (mapcar #'evaluate-fitness (population ea)))

(defmethod select-parents ((ea ea))
  ;; Find the MOST useful.
  (setf (population ea) (sort (population ea) #'> :key 'fitness))
  (setf (parents ea)
        (subseq (population ea)
                0
                (floor (* (length (population ea))
                          (crossover-prob ea))))))

(defmethod recombine-parents ((ea ea))
  (with-slots (children parents) ea
    (setf children nil)
    (loop while (< (length children) (length parents)) do
          (let (father mother child)
            (setf father (random-element parents))
            (setf mother (random-element parents))
            (setf child (crossover father mother))
            (noisy "Father ~A.~%" father)
            (noisy "Mother ~A.~%" mother)
            (noisy "Child: ~A.~%" child)
            (push child children)))))

(defmethod mutate-children ((ea ea))
  (mapcar #'mutate! (children ea)))

(defmethod integrate-children ((ea ea))
  (mapcar (lambda (child)
            (noisy "Child: ~A.~%" child))
          (children ea))
  (setf (population ea) (concatenate 'list (population ea) (children ea))
        (children ea) nil))

(defmethod select-survivors ((ea ea))
  ;; Find the LEAST useful.
  (setf (population ea) (sort (population ea) #'< :key 'fitness))
  (loop while (< (max-pop-size ea) (length (population ea))) do
        (noisy "Eliminating ~A.~%" (first (population ea)))
        (pop (population ea))))

(defmethod terminate-evolution? ((ea ea))
  (< 50 (generation ea)))

(defmethod evolve ((ea ea))
  (noisy "Starting evolution.~%")
  (noisy "Initializing a random population.~%")
  (initialize-random-population ea)
  (noisy "Evaluating the fitness.~%")
  (evaluate-fitness ea)
  (loop until (terminate-evolution? ea) do
        (noisy "Generation ~A loop.~%" (generation ea))
        (noisy "Selecting parents.~%")
        (select-parents ea)
        (mapcar (lambda (parent)
                  (noisy "Parent ~A.~%" parent))
                (parents ea))
        (noisy "Recombining parents.~%")
        (recombine-parents ea)
        (noisy "Mutating children.~%")
        (mutate-children ea)
        (noisy "Integrating the children into the main population.~%")
        (integrate-children ea)
        (noisy "Evaluating the fitness.~%")
        (evaluate-fitness ea)
        (noisy "Selecting survivors.~%")
        (select-survivors ea)
        (setf (population ea) (sort (population ea) #'> :key #'fitness))
        (noisy "Most fit member: ~A.~%" (first (population ea)))
        (incf (generation ea))))

(defmethod fitness-function (individual (rastrigin2d-ea rastrigin2d-ea))
  ;; The Rastrigin function is optimal at 0, minimizing.  We code the rest of
  ;; the EA to assume that a larger fitness value implies a more fit individual,
  ;; so we just invert the sign of the final result of the Rastrigin function to
  ;; produce the fitness.
  (- (rastrigin (coefficient-a (ea individual))
                (mapcar #'value (coerce (genotype individual) 'list)))))

(defmethod fitness-function (individual (easy2d-ea easy2d-ea))
  (let ((v (mapcar #'value (coerce (genotype individual) 'list))))
    (- (* (expt 2.7 (first v))
          (expt 4.2 (second v))))))

(defvar *current-ea* nil)
(defun run-the-ea ()
  (setf *random-state* (make-random-state t))
  (setf *current-ea* (make-instance 'rastrigin2d-ea))
  ;(setf *current-ea* (make-instance 'easy2d-ea))
  (evolve *current-ea*))
