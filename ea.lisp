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

;;; Cf. http://www.necessaryandsufficient.net/2010/12/classical-test-functions-for-genetic-algorithms/ 
(defvar *rastrigin-lower* -5.12)
(defvar *rastrigin-upper* 5.12)
(defun rastrigin (a x)
  (assert (floatp a))
  (assert (listp x))
  (mapcar (lambda (xi)
	    (assert (floatp xi))
	    (assert (< *rastrigin-lower* xi *rastrigin-upper*))))
  (let ((n (length x)))
    (+ (* a n)
       (reduce #'+ (mapcar (lambda (xi)
			     (- (expt xi 2.0)
				(* a (cos (* 2 pi xi)))))
			   x)))))

(defun bounded (min value max)
  (let ((a (min min max))
	(z (max min max)))
    (min (max a value) z)))

(defclass gene ()
  ((value :accessor value :initarg :value
	  :type float
	  :initform nil)
   (lower-bound :accessor lower-bound :initarg :lower-bound
		:type float
		:initform nil)
   (upper-bound :accessor upper-bound :initarg :upper-bound
		:type float
		:initform nil)
   (individual :accessor individual :initarg :individual
	       :type individual
	       :initform nil)))

(defclass individual ()
  ((genotype :accessor genotype :initarg :genotype
	     :type (vector gene)
	     :initform nil)
   (fitness :accessor fitness :initarg :fitness
	    :type float
	    :initform nil)
   (ea :accessor ea :initarg :ea
       :type ea
       :initform nil)))

(defclass ea ()
  ((population :accessor population :initarg :population
	       :type list
	       :initform nil)
   (environ-lower :accessor environ-lower :initarg :environ-lower
		  :type (vector float)
		  :initform nil)
   (environ-upper :accessor environ-upper :initarg :environ-upper
		  :type (vector float)
		  :initform nil)
   (min-pop-size :accessor min-pop-size :initarg :min-pop-size
		 :type (integer 1 *)
		 :initform 1000)
   (max-pop-size :accessor max-pop-size :initarg :max-pop-size
		 :type (integer 1 *)
		 :initform 2000)
   (mutation-prob :accessor mutation-prob :initarg :mutation-prob
		  :type (float 0.0 1.0)
		  :initform 0.1)
   (mutation-factor :accessor mutation-factor :initarg :mutation-factor
		    :type (float 0.0 1.0)
		    :initform 0.05)
   (crossover-prob :accessor crossover-prob :initarg :crossover-prob
		   :type (float 0.0 1.0)
		   :initform 0.1)))

(defgeneric duplicate (thing))
(defgeneric mutate? (individual))
(defgeneric crossover? (individual))
(defgeneric mutate (thing))
(defgeneric mutate! (thing))
(defgeneric crossover (father mother))
(defgeneric make-random-individual (ea))
(defgeneric add-random-individual (ea))

(defmethod duplicate ((gene gene))
  (with-slots (value lower-bound upper-bound individual) gene
    (make-instance 'gene
		   :value value :lower-bound lower-bound :upper-bound upper-bound
		   :individual individual)))

(defmethod duplicate ((individual individual))
  (with-slots (genotype fitness ea) individual
    (make-instance 'individual :genotype genotype :fitness fitness :ea ea)))

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
      (loop for i from 0 to (length genotype) do
	   (mutate! (aref genotype i)))))
  individual)

(defmethod mutate ((individual individual))
  (mutate! (duplicate individual)))

(defmethod crossover ((father gene) (mother gene))
  (let ((result (duplicate gene)))
    (when (= 0 (random 1))
      (setf (value father) (value mother)))
    result))

(defmethod crossover ((father individual) (mother individual))
  (let ((result (duplicate individual)))
    (loop for i from 0 to (length (genotype father)) do
	 (setf (aref (genotype result) i)
	       (mutate (genotype father) (genotype mother))))
    result))

(defmethod make-random-individual ((ea ea))
  (make-instance 'individual
		 :genotype
		 (coerce (loop for i from 0 to (length (environ-lower ea)) collect
			      (let* ((l (aref (environ-lower ea) i))
				     (u (aref (environ-upper ea) i))
				     (a (min l u))
				     (z (max l u)))
				(+ a (random (- z a)))))
			 'vector)
		 :ea ea))

(defmethod add-random-individual ((ea ea))
  (setf (population ea)
	(cons (make-random-individual ea) (population ea))))

(defmethod initialize-random-population ((ea ea))
  (setf (population ea) nil)
  (while (< (length (population ea)) (min-pop-size ea))
    (add-random-individual ea)))
