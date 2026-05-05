(ns kyleerhabor.2d-wavelets-and-image-compression.core
  (:require
   [clojure.math :as math]
   [clojure.string :as str]
   [clojure.test :refer [is]]
   [clojure.tools.cli :as cli]
   [clojure.tools.logging :as log]))

(def sqrt-two (math/sqrt 2))
;; For a grayscale image, 2^8 - 1 = 255
(def default-image-bits-per-pixel 8)
;; 16 * 16 = 256
(def default-image-size 16)

;; s[k, j] = (s[k + 1, 2j] + s[k + 1, 2j + 1]) / √2
(defn approximation [s k j]
  (let [jj (* 2 j)
        x (/ (+ (get s jj) (get s (inc jj))) sqrt-two)]
    (log/trace (str "approximation: k = " k ", j = " j ", s[k, j] = " x))
    x))

;; a[k, j] = (s[k + 1, 2j] - s[k + 1, 2j + 1]) / √2
(defn detail [s k j]
  (let [jj (* 2 j)
        x (/ (- (get s jj) (get s (inc jj))) sqrt-two)]
    (log/trace (str "detail: k = " k ", j = " j ", a[k, j] = " x))
    x))

;; s[k + 1, 2j] = (s[k, j] + a[k, j]) / √2
(defn invert-1 [s a k j]
  (let [x (/ (+ s a) sqrt-two)]
    (log/trace (str "invert: k = " k ", j = " j ", s[k + 1, 2j] = " x))
    x))

;; s[k + 1, 2j + 1] = (s[k, j] - a[k, j]) / √2
(defn invert-2 [s a k j]
  (let [x (/ (- s a) sqrt-two)]
    (log/trace (str "invert: k = " k ", j = " j ", s[k + 1, 2j] = " x))
    x))

(defn power-of-two? [x]
  ;; https://stackoverflow.com/a/600306
  (and (pos? x) (zero? (bit-and x (dec x)))))

(defn digit-count [n]
  ;; https://stackoverflow.com/a/40322349
  (if (zero? n)
    1
    (int (inc (math/log10 n)))))

(defn image-row-name [row ds]
  (let [precision 2
        ;; The 2 is for the minus sign and decimal point.
        padding (+ 2 ds precision)]
    (str/join " " (map #(format (str "%" padding "." precision "f") %) row))))

(defn image-name [{:keys [data size]}]
  (let [ds (digit-count (int (reduce max 0 data)))
        rows (->> data 
               (partition size)
               (map #(image-row-name % ds))
               (str/join "\n "))
        result (str "[" rows "]")]
    result))

(defn transform-row [row size k]
  (let [js (range (quot size 2))
        s (map #(approximation row k %) js)
        a (map #(detail row k %) js)
        result (concat s a)]
    result))

(defn transform-rows [{:keys [data size]} k]
  (let [data (->> data
               (partition size)
               (map #(transform-row (vec %) size k))
               (reduce concat))
        image {:data data
               :size size}]
    image))

(defn transpose-row [v size i] 
  (->> i
    (iterate #(+ % size))
    (take size)
    (map #(get v %))))

(defn transpose [{:keys [data size]}]
  (let [v (vec data)
        data (->> (range size)
               (map #(transpose-row v size %))
               (reduce concat))
        image {:data data
               :size size}]
    image))

(defn transform-columns [image k]
  (transpose (transform-rows (transpose image) k)))

(def percentile-90-threshold 0.9)

(defn percentile-90 [s] 
  (let [ss (sort (map abs s))
        i (int (math/floor (* percentile-90-threshold (count s))))
        threshold (nth ss i)]
    threshold))

(defn edge-axes-2 [{:keys [l h]} [low high]]
  {:l (concat l low)
   :h (concat h high)})

(defn edge-axes [half h]
  (->> h
    (map #(split-at half %))
    (reduce edge-axes-2 {:l [] :h []})))

(defn edges [{:keys [data size]}]
  (let [half (/ size 2)
        rows (partition size data)
        [l h] (split-at half rows)
        {ll :l lh :h} (edge-axes half l)
        {hl :l hh :h} (edge-axes half h)
        result {:ll ll
                :lh lh
                :hl hl
                :hh hh}]
    result))

(defn edges->image [{:keys [size]} {:keys [ll lh hl hh]}]
  (let [half (/ size 2)
        l (concat ll hl)
        h (concat lh hh)
        data (reduce concat (interleave (partition half l) (partition half h)))
        image {:data data
               :size size}]
    image))

(defn apply-threshold [threshold x]
  (if (< (abs x) threshold)
    0.0
    x))

(defn invert-transform-row-2 [[s a] k j]
  {:s (invert-1 s a k j)
   :a (invert-2 s a k j)})

(defn invert-transform-row-2-2 [result x]
  {:s (conj (:s result) (:s x))
   :a (conj (:a result) (:a x))})

(defn invert-transform-row [row k]
  (let [[s a] (split-at (/ (count row) 2) row)
        {:keys [s a]} (->> (interleave s a)
                        (partition 2)
                        (map-indexed #(invert-transform-row-2 %2 k %1))
                        (reduce invert-transform-row-2-2 {:s []
                                                          :a []}))
        result (interleave s a)]
    result))

(defn invert-transform-rows [{:keys [data size]} k]
  (let [data (->> data
               (partition size)
               (map #(invert-transform-row % k))
               (reduce concat))
        image {:data data
               :size size}]
    image))

(defn invert-transform-columns [image k]
  (transpose (invert-transform-rows (transpose image) k)))

(defn rmse [image image']
  (let [data (:data image)
        data' (:data image')
        n (count data)
        sum-squared-diff (reduce + (map #(math/pow (- %1 %2) 2) data data'))]
    (math/sqrt (/ sum-squared-diff n))))

(defn theorem [{:keys [size]
                :as image}]
  (println (image-name image))
  (let [k (dec (Long/numberOfTrailingZeros size))
        image' (transform-rows image k)
        _ (println (image-name image'))
        image' (transform-columns image' k)
        _ (println (image-name image'))
        {:keys [ll lh hl hh]} (edges image')
        threshold (percentile-90 (concat lh hl hh))
        lh (map #(apply-threshold threshold %) lh)
        hl (map #(apply-threshold threshold %) hl)
        hh (map #(apply-threshold threshold %) hh)
        image' (edges->image image' {:ll ll :lh lh :hl hl :hh hh})
        _ (println (image-name image'))
        image' (invert-transform-columns image' k)
        _ (println (image-name image'))
        image' (invert-transform-rows image' k)
        _ (println (image-name image'))]
    (print "RMSE:" (format "%.2f" (rmse image image')))
    image'))

(defn ->image [bits-per-pixel size]
  {:data (repeatedly (* size size) #(double (rand-int (math/pow 2 bits-per-pixel))))
   :size size})

(comment
  (defn test-image [image data]
    (is (= (map math/round (:data image)) data)))
  
  (def image {:data [194.0 157.0 221.0 164.0
                      13.0 209.0 182.0  88.0
                     152.0 200.0  14.0 139.0
                       0.0 199.0   8.0  48.0]
              :size 4})
  (test-image (theorem image) [201  85 164 164
                                85 201 164 164
                                76 199  52  52
                                76 199  52  52])
  
  (def image {:data [146.0  70.0 159.0  44.0   3.0  28.0 234.0 238.0
                      24.0 101.0  30.0  87.0 247.0 209.0 194.0  42.0
                      72.0 186.0 111.0  62.0 165.0 117.0 152.0 113.0
                     115.0  58.0 176.0  83.0 182.0 129.0   3.0 155.0
                      25.0 201.0 247.0  43.0   7.0  84.0  79.0 213.0
                     248.0 156.0  13.0 169.0  31.0 245.0 252.0  34.0
                     161.0   6.0 122.0  66.0   8.0 118.0  37.0  32.0
                      60.0 171.0  31.0  83.0  77.0 241.0  54.0  49.0]
              :size 8})
  (test-image (theorem image) [ 85  85  80  80  15  15 177 177
                                85  85  80  80 228 228 177 177
                               108 108 108 108 148 148 106 106
                               108 108 108 108 148 148 106 106
                               157 157 208  28  19 164  57 232
                               157 157  28 208  19 164 232  57
                                99  99  75  75  42 179  43  43
                                99  99  75  75  42 179  43  43])
  
  (theorem (->image default-image-bits-per-pixel default-image-size))
  
  ;; [a b c d
  ;;  e f g h
  ;;  i j k l
  ;;  m n o p]
  ;; 
  ;; [a' b* c' d*
  ;;  e' f* g' h*
  ;;  i' j* k' l*
  ;;  m' n* o' p*]
  ;; 
  ;; [a'' b*' c'' d*'
  ;;  e'* f** g'* h**
  ;;  i'' j*' k'' l*'
  ;;  m'* n** o'* p**]
  )

(def cli-options
  [["-h" "--help"]
   [nil "--bits-per-pixel NUMBER" nil
    :default default-image-bits-per-pixel
    :parse-fn parse-long
    :validate [#(not (neg? %)) "Must not be negative"]]
   [nil "--size NUMBER" nil
    :default default-image-size
    :parse-fn parse-long
    :validate [power-of-two? "Must be a power of 2"]]])

(defn -main [& args]
  (let [{:keys [summary errors]
         {:keys [help bits-per-pixel size]} :options} (cli/parse-opts args cli-options)]
    (cond
      help (println summary)
      errors (println (str/join "\n" errors))
      :else (theorem (->image bits-per-pixel size)))))