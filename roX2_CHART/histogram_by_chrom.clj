; Display the distribution of reads per chromosome from an aligned BAM file.
;
; Usage:
;   cljr run histogram_by_chrom.clj <in BAM file> [<chromosomes to remove>]
; 
; Requires:
;   Picard
;   Incanter

(import '[net.sf.samtools SAMFileReader])
(use '[clojure.java.io]
     '[clojure.contrib.command-line])
(require '[clojure.string :as s])
(require '[incanter core charts])

(defn bam-chromosomes [bam-iter]
  "Lazy sequence of aligned chromosome names from a BAM file."
  (lazy-seq
    (if (.hasNext bam-iter)
      (let [rec (.next bam-iter)]
        (cons 
          (if-not (.getReadUnmappedFlag rec) (.getReferenceName rec))
          (bam-chromosomes bam-iter))))))

(defn chrom-count [chrom-coll]
  "Provide a map of counts for each chromosome in a collection."
  (reduce (fn [counts c] (assoc counts c (inc (get counts c 0))))
          {} (remove nil? chrom-coll)))

(defn sizes-from-bam [bam-reader]
  "Provide a map of sizes for each chromosome from the BAM file header."
  (reduce (fn [sizes [n l]] (assoc sizes n l))
          {} (for [rec (-> bam-reader 
                         .getFileHeader .getSequenceDictionary .getSequences)]
              [(.getSequenceName rec) (.getSequenceLength rec)])))

(defn bam-to-histogram [in-file to-remove]
  "Produce a histogram of counts per chromosome from an aligned BAM file."
  (let [bam-reader (SAMFileReader. (file in-file))
        bam-iter (.iterator bam-reader)
	base-name (last (s/split in-file #"/"))
        out-file (s/join "" [(first (s/split base-name #"\.\w+$")) "-chrom.png"])
        size-map (sizes-from-bam bam-reader)
        count-map (chrom-count (bam-chromosomes bam-iter))
        chroms (sort (remove #(contains? to-remove %) (keys count-map)))
        mb-scale #(* 1e6 %)
        reads (for [c chroms] (/ (mb-scale (get count-map c)) (get size-map c)))]
    (doto (incanter.charts/bar-chart chroms reads
                                     :vertical false
                                     :title base-name
                                     :x-label ""
                                     :y-label "reads per Mb")
      (incanter.core/save out-file :width 500 :height 400))))

(when *command-line-args*
  (let [in-file (first *command-line-args*)
        to-remove (if-not (nil? (second *command-line-args*))
                    (s/split (second *command-line-args*) #",")
                    [])]
    (bam-to-histogram in-file (set to-remove))))
