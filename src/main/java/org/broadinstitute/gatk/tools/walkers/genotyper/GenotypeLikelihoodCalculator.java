/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system ("PHONE-HOME") which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE'S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2016 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gatk.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.collections.IntMaxHeap;
import org.broadinstitute.gatk.utils.genotyper.ReadLikelihoods;

/**
 * Helper to calculate genotype likelihoods given a ploidy and an allele count (number of possible distinct alleles).
 *
 * <p>
 *     Notice that for performance this class is thread-unsafe an so it cannot be shared between thread in a multi-thread run.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class GenotypeLikelihoodCalculator {

    /**
     * Maximum number of components (or distinct alleles) for any genotype with this calculator ploidy and allele count.
     */
    private int maximumDistinctAllelesInGenotype;

    /**
     * Offset table for this calculator.
     *
     * <p>
     *     This is a shallow copy of {@link GenotypeLikelihoodCalculators#alleleFirstGenotypeOffsetByPloidy} when the calculator was created
     *     thus it follows the same format as that array. Please refer to its documentation.
     * </p>
     *
     * <p>You can assume that this offset table contain at least (probably more) the numbers corresponding to the allele count and ploidy for this calculator.
     * However since it might have more than that and so you must use {@link #alleleCount} and {@link #ploidy} when
     * iterating through this array rather that its length or the length of its components.</p>.
     */
    private final int[][] alleleFirstGenotypeOffsetByPloidy;

    /**
     * Genotype table for this calculator.
     *
     * <p>It is ensure that it contains all the genotypes for this calculator ploidy and allele count, maybe more. For
     * that reason you must use {@link #genotypeCount} when iterating through this array and not relay on its length.</p>
     */
    private final GenotypeAlleleCounts[] genotypeAlleleCounts;

    /**
     * Number of genotypes given this calculator {@link #ploidy} and {@link #alleleCount}.
     */
    private final int genotypeCount;

    /**
     * Number of genotyping alleles for this calculator.
     */
    private final int alleleCount;

    /**
     * Ploidy for this calculator.
     */
    private final int ploidy;

    /**
     * Max-heap for integers used for this calculator internally.
     */
    private final IntMaxHeap alleleHeap;

    /**
     * Cache of the last genotype-allele-count requested using {@link #genotypeAlleleCountsAt(int)}, when it
     * goes beyond the maximum genotype-allele-count static capacity. Check on that method documentation for details.
     */
    private transient GenotypeAlleleCounts lastOverheadCounts;

    /**
     * Buffer used as a temporary container for likelihood components for genotypes stratified by alleles, allele frequency and reads.
     *
     * <p>To improve performance we use a 1-dimensional array to implement a 3-dimensional one as some of those dimension
     * have typically very low depths (allele and allele frequency)</p>
     *
     * <p>
     *     The value contained in position <code>[a][f][r] == log10Lk(read[r] | allele[a]) + log10(f) </code>. Exception is
     *     for f == 0 whose value is undefined (in practice 0.0) and never used.
     * </p>
     *
     * <p>
     *     It is indexed by read, then by allele and then by the number of copies of the allele. For the latter
     *     there are as many entries as the ploidy of the calculator + 1 (to accommodate zero copies although is
     *     never used in practice).
     * </p>
     */
    private double[] readAlleleLikelihoodByAlleleCount = null;

    /**
     * Buffer used as a temporary container for likelihood components for genotypes stratified by reads.
     *
     * <p>
     *     It is indexed by genotype index and then by read index. The read capacity is increased as needed by calling
     *     {@link #ensureReadCapacity(int) ensureReadCapacity}.
     * </p>
     */
    private final double[][] readLikelihoodsByGenotypeIndex;

    /**
     * Indicates how many reads the calculator supports.
     *
     * <p>This figure is increased dynamically as per the
     * calculation request calling {@link #ensureReadCapacity(int) ensureReadCapacity}.<p/>
     */
    private int readCapacity = -1;

    /**
     * Caches the log10 of the first few integers up to the ploidy supported by the calculator.
     * <p>This is in fact a shallow copy if {@link GenotypeLikelihoodCalculators#ploidyLog10}</p> and is not meant to be modified by
     * this class. </p>
     */
    private final double[] ploidyLog10;

    /**
     * Buffer field use as a temporal container for sorted allele counts when calculating the likelihood of a
     * read in a genotype.
     * <p>
     *      This array follows the same format as {@link GenotypeAlleleCounts#sortedAlleleCounts}. Each component in the
     *      genotype takes up two positions in the array where the first indicate the allele index and the second its frequency in the
     *      genotype. Only non-zero frequency alleles are represented, taking up the first positions of the array.
     * </p>
     *
     * <p>
     *     This array is sized so that it can accommodate the maximum possible number of distinct alleles in any
     *     genotype supported by the calculator, value stored in {@link #maximumDistinctAllelesInGenotype}.
     * </p>
     */
    private final int[] genotypeAllelesAndCounts;

    /**
     * Buffer field use as a temporal container for component likelihoods when calculating the likelihood of a
     * read in a genotype. It is stratified by read and the allele component of the genotype likelihood... that is
     * the part of the likelihood sum that correspond to a particular allele in the genotype.
     *
     * <p>
     *     It is implemented in a 1-dimensional array since typically one of the dimensions is rather small. Its size
     *     is equal to {@link #readCapacity} times {@link #maximumDistinctAllelesInGenotype}.
     * </p>
     *
     * <p>
     *     More concretely [r][i] == log10Lk(read[r] | allele[i]) + log(freq[i]) where allele[i] is the ith allele
     *     in the genotype of interest and freq[i] is the number of times it occurs in that genotype.
     * </p>
     */
    private double[] readGenotypeLikelihoodComponents;

    /**
     * Creates a new calculator providing its ploidy and number of genotyping alleles.
     */
    protected GenotypeLikelihoodCalculator(final int ploidy, final int alleleCount,
                                           final int[][] alleleFirstGenotypeOffsetByPloidy,
                                           final GenotypeAlleleCounts[][] genotypeTableByPloidy,
                                           final double[] ploidyLog10) {

        this.alleleFirstGenotypeOffsetByPloidy = alleleFirstGenotypeOffsetByPloidy;
        genotypeAlleleCounts = genotypeTableByPloidy[ploidy];
        this.alleleCount = alleleCount;
        this.ploidy = ploidy;
        genotypeCount = this.alleleFirstGenotypeOffsetByPloidy[ploidy][alleleCount];

        alleleHeap = new IntMaxHeap(ploidy);
        readLikelihoodsByGenotypeIndex = new double[genotypeCount][];
        this.ploidyLog10 = ploidyLog10;
        // The number of possible components is limited by distinct allele count and ploidy.
        maximumDistinctAllelesInGenotype = Math.min(ploidy, alleleCount);
        genotypeAllelesAndCounts = new int[maximumDistinctAllelesInGenotype*2];
    }

    /**
     * Makes sure that temporal arrays and matrices are prepared for a number of reads to process.
     * @param requestedCapacity number of read that need to be processed.
     */
    public void ensureReadCapacity(final int requestedCapacity) {
        if (requestedCapacity < 0)
            throw new IllegalArgumentException("illegal capacity value");
        if (readCapacity == -1) { // first time call.
            final int minimumCapacity = Math.max(requestedCapacity,10); // Never go too small, 10 is the minimum.
            readAlleleLikelihoodByAlleleCount = new double[minimumCapacity * alleleCount * (ploidy+1)];
            for (int i = 0; i < genotypeCount; i++)
                readLikelihoodsByGenotypeIndex[i] = new double[minimumCapacity];
            readGenotypeLikelihoodComponents = new double[ploidy * minimumCapacity];
            readCapacity = minimumCapacity;
        } else if (readCapacity < requestedCapacity) {
            final int doubleCapacity = (requestedCapacity << 1);
            readAlleleLikelihoodByAlleleCount = new double[doubleCapacity * alleleCount * (ploidy+1)];
            for (int i = 0; i < genotypeCount; i++)
                readLikelihoodsByGenotypeIndex[i] = new double[doubleCapacity];
            readGenotypeLikelihoodComponents = new double[maximumDistinctAllelesInGenotype * doubleCapacity];
            readCapacity = doubleCapacity;
        }
    }

    /**
     * Give a list of alleles, returns the likelihood array index.
     *
     * <p>This operation is <b>thread-unsafe</b>.</p>
     *
     * @param alleleIndices the indices of the alleles in the genotype, there should be as many repetition of an
     *                      index as copies of that allele in the genotype. Allele indices do not need to be sorted in
     *                      any particular way.
     *
     * @return never {@code null}.
     */
    public int allelesToIndex(final int... alleleIndices) {
        // Special case ploidy == 0.
        if (ploidy == 0) return 0;

        alleleHeap.clear();
        alleleHeap.add(alleleIndices);
        return alleleHeapToIndex();
    }

    /**
     * Returns the number of possible genotypes given ploidy and the maximum allele index.
     * @return never {@code null}.
     */
    public int genotypeCount()  {
        return genotypeCount;
    }

    /**
     * Returns the genotype associated to a particular likelihood index.
     *
     * <p>If {@code index} is larger than {@link GenotypeLikelihoodCalculators#MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY},
     *  this method will reconstruct that genotype-allele-count iteratively from the largest strongly referenced count available.
     *  or the last requested index genotype.
     *  </p>
     *
     * <p> Therefore if you are iterating through all genotype-allele-counts you should do sequentially and incrementally, to
     * avoid a large efficiency drop </p>.
     *
     * @param index query likelihood-index.
     * @return never {@code null}.
     */
    public GenotypeAlleleCounts genotypeAlleleCountsAt(final int index) {
        if (index < 0 || index >= genotypeCount)
            throw new IllegalArgumentException("invalid likelihood index: " + index + " >= " + genotypeCount
                    + " (genotype count for nalleles = " + alleleCount + " and ploidy " + ploidy );
        if (index < GenotypeLikelihoodCalculators.MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY)
            return genotypeAlleleCounts[index];
        else if (lastOverheadCounts == null || lastOverheadCounts.index() > index) {
            final GenotypeAlleleCounts result = genotypeAlleleCounts[GenotypeLikelihoodCalculators.MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY - 1].clone();
            result.increase(index - GenotypeLikelihoodCalculators.MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY + 1);
            lastOverheadCounts = result;
            return result.clone();
        } else {
            lastOverheadCounts.increase(index - lastOverheadCounts.index());
            return lastOverheadCounts.clone();
        }
    }

    /**
     * Calculate the likelihoods given the list of alleles and the likelihood map.
     *
     * <p>This operation is <b>thread-unsafe</b>.</p>
     *
     * @param likelihoods the likelihood matrix all alleles vs all reads.
     *
     * @throws IllegalArgumentException if {@code alleleList} is {@code null} or {@code likelihoods} is {@code null}
     *     or the alleleList size does not match the allele-count of this calculator, or there are missing allele vs
     *     read combinations in {@code likelihoods}.
     *
     * @return never {@code null}.
     */
    public <A extends Allele> GenotypeLikelihoods genotypeLikelihoods(final ReadLikelihoods.Matrix<A> likelihoods) {
        if (likelihoods == null)
            throw new IllegalArgumentException("the likelihood map cannot be null");

        if (likelihoods.alleleCount() != alleleCount)
            throw new IllegalArgumentException("mismatch between allele list and alleleCount");


        final int readCount = likelihoods.readCount();


        ensureReadCapacity(readCount);

        /// [x][y][z] = z * LnLk(Read_x | Allele_y)
        final double[] readLikelihoodComponentsByAlleleCount
                = readLikelihoodComponentsByAlleleCount(likelihoods);
        final double[][] genotypeLikelihoodByRead = genotypeLikelihoodByRead(readLikelihoodComponentsByAlleleCount,readCount);

        final double[] readLikelihoodsByGenotypeIndex = genotypeLikelihoods(genotypeLikelihoodByRead, readCount);
        return GenotypeLikelihoods.fromLog10Likelihoods(readLikelihoodsByGenotypeIndex);
    }

    /**
     * Calculates the final genotype likelihood array out of the likelihoods for each genotype per read.
     *
     * @param readLikelihoodsByGenotypeIndex <i>[g][r]</i> likelihoods for each genotype <i>g</i> and <i>r</i>.
     * @param readCount number of reads in the input likelihood arrays in {@code genotypeLikelihoodByRead}.
     * @return never {@code null}, one position per genotype where the <i>i</i> entry is the likelihood of the ith
     *   genotype (0-based).
     */
    public double[] genotypeLikelihoods(final double[][] readLikelihoodsByGenotypeIndex, final int readCount) {
        final double[] result = new double[genotypeCount];
        final double denominator = readCount * ploidyLog10[ploidy]; // instead of dividing each read likelihood by ploidy
         // ( so subtract log10(ploidy) )  we multiply them all and the divide by ploidy^readCount (so substract readCount * log10(ploidy) )
        for (int g = 0; g < genotypeCount; g++) {
            final double[] likelihoodsByRead = readLikelihoodsByGenotypeIndex[g];
            //double s = - denominator;
            double s = 0.0;
            System.out.println("denominator = " + denominator);
            for (int r = 0; r < readCount; r++) {
                s += likelihoodsByRead[r];
                System.out.println("genotypeLikelihoods:" + r  + ":" + likelihoodsByRead[r]);
            }
            //result[g] = s;
            result[g] = s - denominator;
            System.out.println("genotypeLikelihoods(output):" + g + "," + s + "," + denominator);
        }
        return result;
    }

    /**
     * Calculates the likelihood component of each read on each genotype.
     *
     * @param readLikelihoodComponentsByAlleleCount [a][f][r] likelihood stratified by allele <i>a</i>, frequency in genotype <i>f</i> and
     *                                              read <i>r</i>.
     * @param readCount number of reads in {@code readLikelihoodComponentsByAlleleCount}.
     * @return never {@code null}.
     */
    public double[][] genotypeLikelihoodByRead(final double[] readLikelihoodComponentsByAlleleCount, final int readCount) {

        // Here we don't use the convenience of {@link #genotypeAlleleCountsAt(int)} within the loop to spare instantiations of
        // GenotypeAlleleCounts class when we are dealing with many genotypes.
        GenotypeAlleleCounts alleleCounts = genotypeAlleleCounts[0];

        for (int genotypeIndex = 0; genotypeIndex < genotypeCount; genotypeIndex++) {
            final double[] readLikelihoods = this.readLikelihoodsByGenotypeIndex[genotypeIndex];
            final int componentCount = alleleCounts.distinctAlleleCount();
            System.out.println("component count" + componentCount + ", genotypeIndex= " + genotypeIndex);
            switch (componentCount) {
                case 1: //
                    singleComponentGenotypeLikelihoodByRead(alleleCounts, readLikelihoods, readLikelihoodComponentsByAlleleCount, readCount);
                    break;
                case 2:
                    twoComponentGenotypeLikelihoodByRead(alleleCounts,readLikelihoods,readLikelihoodComponentsByAlleleCount, readCount);
                    break;
                default:
                    manyComponentGenotypeLikelihoodByRead(alleleCounts,readLikelihoods,readLikelihoodComponentsByAlleleCount, readCount);
            }
            if (genotypeIndex < genotypeCount - 1)
                alleleCounts = nextGenotypeAlleleCounts(alleleCounts);
        }
        return readLikelihoodsByGenotypeIndex;
    }

    private GenotypeAlleleCounts nextGenotypeAlleleCounts(final GenotypeAlleleCounts alleleCounts) {
        final int index = alleleCounts.index();
        final GenotypeAlleleCounts result;
        final int cmp = index - GenotypeLikelihoodCalculators.MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY + 1;
        if (cmp < 0)
            result = genotypeAlleleCounts[index + 1];
        else if (cmp == 0) {
            result = genotypeAlleleCounts[index].clone();
            result.increase();
        } else {
            alleleCounts.increase();
            result = alleleCounts;
        }
        return result;
    }

    /**
     * General genotype likelihood component by thread calculator. It does not make any assumption in the exact
     * number of alleles present in the genotype.
     */
    private void manyComponentGenotypeLikelihoodByRead(final GenotypeAlleleCounts genotypeAlleleCounts,
                                                       final double[] likelihoodByRead,
                                                       final double[]readLikelihoodComponentsByAlleleCount,
                                                       final int readCount) {

        // First we collect the allele likelihood component for all reads and place it
        // in readGenotypeLikelihoodComponents for the final calculation per read.
        System.out.println("manyComponentGenotypeLikelihoodByRead");
        genotypeAlleleCounts.copyAlleleCounts(genotypeAllelesAndCounts,0);
        final int componentCount = genotypeAlleleCounts.distinctAlleleCount();
        final int alleleDataSize = (ploidy + 1) * readCount;
        for (int c = 0,cc = 0; c < componentCount; c++) {
            final int alleleIndex = genotypeAllelesAndCounts[cc++];
            final int alleleCount = genotypeAllelesAndCounts[cc++];
            // alleleDataOffset will point to the index of the first read likelihood for that allele and allele count.
            int alleleDataOffset = alleleDataSize * alleleIndex + alleleCount * readCount;
            for (int r = 0, readDataOffset = c; r < readCount; r++, readDataOffset += maximumDistinctAllelesInGenotype)
                readGenotypeLikelihoodComponents[readDataOffset] = readLikelihoodComponentsByAlleleCount[alleleDataOffset++];
        }

        // Calculate the likelihood per read.
        for (int r = 0, readDataOffset = 0; r < readCount; r++, readDataOffset += maximumDistinctAllelesInGenotype)
            likelihoodByRead[r] = MathUtils.approximateLog10SumLog10(readGenotypeLikelihoodComponents, readDataOffset, readDataOffset + componentCount);
    }

    /**
     * Calculates the likelihood component by read for a given genotype allele count assuming that there are
     * exactly two alleles present in the genotype (with arbitrary non-zero counts each).
     */
    private void twoComponentGenotypeLikelihoodByRead(final GenotypeAlleleCounts genotypeAlleleCounts,
                                                      final double[] likelihoodByRead,
                                                      final double[] readLikelihoodComponentsByAlleleCount,
                                                      final int readCount) {

        final int allele0 = genotypeAlleleCounts.alleleIndexAt(0);
        final int freq0 = genotypeAlleleCounts.alleleCountAt(0);
        final int allele1 = genotypeAlleleCounts.alleleIndexAt(1);
        final int freq1 = ploidy - freq0; // no need to get it from genotypeAlleleCounts.
        int allele0LnLkOffset = readCount * ((ploidy + 1) * allele0 + freq0);
        int allele1LnLkOffset = readCount * ((ploidy + 1) * allele1 + freq1);
        System.out.println("twoComponentGenotypeLikelihoodByRead " + readCount + "," + allele0 + "," +
                freq0 + "," + allele1 + "," + freq1 + "," + allele0LnLkOffset + "," + allele1LnLkOffset + ", ploidy=" + ploidy);
        for (int r = 0; r < readCount; r++) {
            System.out.println("twoComponentGenotypeLikelihoodByRead(cycle):" + allele0LnLkOffset +
                    ":" + allele1LnLkOffset + ":" + readLikelihoodComponentsByAlleleCount[allele0LnLkOffset]
                    + ":" + readLikelihoodComponentsByAlleleCount[allele1LnLkOffset] + " = " +
                    MathUtils.approximateLog10SumLog10(readLikelihoodComponentsByAlleleCount[allele0LnLkOffset],readLikelihoodComponentsByAlleleCount[allele1LnLkOffset]));
            final double lnLk0 = readLikelihoodComponentsByAlleleCount[allele0LnLkOffset++];
            final double lnLk1 = readLikelihoodComponentsByAlleleCount[allele1LnLkOffset++];
            likelihoodByRead[r] = MathUtils.approximateLog10SumLog10(lnLk0,lnLk1);
        }
    }

    /**
     * Calculates the likelihood component by read for a given genotype allele count assuming that there are
     * exactly one allele present in the genotype.
     */
    private void singleComponentGenotypeLikelihoodByRead(final GenotypeAlleleCounts genotypeAlleleCounts,
                                                         final double[] likelihoodByRead,
                                                         final double[] readLikelihoodComponentsByAlleleCount,
                                                         final int readCount) {
        final int allele = genotypeAlleleCounts.alleleIndexAt(0);
        System.out.println("singleComponentGenotypeLikelihoodByRead:" + allele + ":" + ploidy + ":" + genotypeAlleleCounts + ":" + readCount);

        // the count of the only component must be = ploidy.
        int offset = (allele * (ploidy + 1) + ploidy) * readCount;
        for (int r = 0; r < readCount; r++) {
            System.out.println("singleComponentGenotypeLikelihoodByRead" + r + ":" + offset + ":" + readLikelihoodComponentsByAlleleCount[offset]);
            likelihoodByRead[r] =
                    readLikelihoodComponentsByAlleleCount[offset++];
        }
    }

    /**
     * Returns a 3rd matrix with the likelihood components.
     *
     * <pre>
     *     result[y][z][x] :=  z * lnLk ( read_x | allele_y ).
     * </pre>
     *
     * @return never {@code null}.
     */
    public <A extends Allele> double[] readLikelihoodComponentsByAlleleCount(final ReadLikelihoods.Matrix<A> likelihoods) {
        final int readCount = likelihoods.readCount();
        final int alleleDataSize = readCount * (ploidy + 1);

        // frequency1Offset = readCount to skip the useless frequency == 0. So now we are at the start frequency == 1
        // frequency1Offset += alleleDataSize to skip to the next allele index data location (+ readCount) at each iteration.
        for (int a = 0, frequency1Offset = readCount; a < alleleCount; a++, frequency1Offset += alleleDataSize) {
            likelihoods.copyAlleleLikelihoods(a, readAlleleLikelihoodByAlleleCount, frequency1Offset);

            // p = 2 because the frequency == 1 we already have it.
            for (int frequency = 2, destinationOffset = frequency1Offset + readCount; frequency <= ploidy; frequency++) {
                final double log10frequency = ploidyLog10[frequency];
                for (int r = 0, sourceOffset = frequency1Offset; r < readCount; r++) {
                    System.out.println("readLikelihoodComponentsByAlleleCount(cycle):" + r + ":" + frequency + ":" +
                            frequency1Offset + ":" + destinationOffset + ":" + sourceOffset +
                    ":" + readAlleleLikelihoodByAlleleCount[destinationOffset] + ":" + readAlleleLikelihoodByAlleleCount[sourceOffset] + ":" + log10frequency);
                    readAlleleLikelihoodByAlleleCount[destinationOffset++] =
                            readAlleleLikelihoodByAlleleCount[sourceOffset++] + log10frequency;
                }
            }
        }
        return readAlleleLikelihoodByAlleleCount;
    }

    /**
     * Returns the ploidy for this genotype likelihood calculator.
     * @return 0 or greater.
     */
    public int ploidy() {
        return ploidy;
    }

    /**
     * Returns the total number of alleles for this genotype calculator.
     * @return the number of alleles considered by this calculator.
     */
    public int alleleCount() {
        return alleleCount;
    }

    /**
     * Returns the likelihood index given the allele counts.
     *
     * @param alleleCountArray the query allele counts. This must follow the format returned by
     *  {@link GenotypeAlleleCounts#copyAlleleCounts} with 0 offset.
     *
     * @throws IllegalArgumentException if {@code alleleCountArray} is not a valid {@code allele count array}:
     *  <ul>
     *      <li>is {@code null},</li>
     *      <li>or its length is not even,</li>
     *      <li>or it contains any negatives,
     *      <li>or the count sum does not match the calculator ploidy,</li>
     *      <li>or any of the alleles therein is negative or greater than the maximum allele index.</li>
     *  </ul>
     *
     * @return 0 or greater but less than {@link #genotypeCount}.
     */
    public int alleleCountsToIndex(final int ... alleleCountArray) {
        if (alleleCountArray == null)
            throw new IllegalArgumentException("the allele counts cannot be null");
        if ((alleleCountArray.length & 1) != 0)
            throw new IllegalArgumentException("the allele counts array cannot have odd length");
        alleleHeap.clear();
        for (int i = 0; i < alleleCountArray.length; i += 2) {
            final int index = alleleCountArray[i];
            final int count = alleleCountArray[i+1];
            if (count < 0)
                throw new IllegalArgumentException("no allele count can be less than 0");
            for (int j = 0; j < count; j++)
                alleleHeap.add(index);
        }
        return alleleHeapToIndex();
    }

    /**
     * Transforms the content of the heap into an index.
     *
     * <p>
     *     The heap contents are flushed as a result, so is left ready for another use.
     * </p>
     *
     * @return a valid likelihood index.
     */
    private int alleleHeapToIndex() {
        if (alleleHeap.size() != ploidy)
            throw new IllegalArgumentException("the sum of allele counts must be equal to the ploidy of the calculator");
        if (alleleHeap.peek() >= alleleCount)
            throw new IllegalArgumentException("invalid allele " + alleleHeap.peek() + " more than the maximum " + (alleleCount - 1));
        int result = 0;
        for (int p = ploidy; p > 0; p--) {
            final int allele = alleleHeap.remove();
            if (allele < 0)
                throw new IllegalArgumentException("invalid allele " + allele + " must be equal or greater than 0 ");
            result += alleleFirstGenotypeOffsetByPloidy[p][allele];
        }
        return result;
    }

    /**
     * Composes a genotype index map given a allele index recoding.
     *
     * @param oldToNewAlleleIndexMap allele recoding. The ith entry indicates the index of the allele in original encoding
     *                               that corresponds to the ith allele index in the final encoding.
     *
     * @throws IllegalArgumentException if this calculator cannot handle the recoding provided. This is
     * the case when either {@code oldToNewAlleleIndexMap}'s length or any of its element (+ 1 as they are 0-based) is larger
     * this calculator's {@link #alleleCount()}. Also if any {@code oldToNewAllelesIndexMap} element is negative.
     *
     * @return never {@code null}.
     */
    public int[] genotypeIndexMap(final int[] oldToNewAlleleIndexMap) {
        if (oldToNewAlleleIndexMap == null)
            throw new IllegalArgumentException("the input encoding array cannot be null");

        final int resultAlleleCount = oldToNewAlleleIndexMap.length;
        if (resultAlleleCount > alleleCount)
            throw new IllegalArgumentException("this calculator does not have enough capacity for handling "
                    + resultAlleleCount + " alleles ");
        final int resultLength = resultAlleleCount == alleleCount
                ? genotypeCount : GenotypeLikelihoodCalculators.genotypeCount(ploidy,resultAlleleCount);

        final int[] result = new int[resultLength];
        final int[] sortedAlleleCounts = new int[Math.max(ploidy,alleleCount) << 1];
        alleleHeap.clear();
        GenotypeAlleleCounts alleleCounts = genotypeAlleleCounts[0];
        for (int i = 0; i < resultLength; i++) {
            genotypeIndexMapPerGenotypeIndex(i,alleleCounts, oldToNewAlleleIndexMap, result, sortedAlleleCounts);
            if (i < resultLength - 1)
                  alleleCounts = nextGenotypeAlleleCounts(alleleCounts);
        }
        return result;
    }

    /**
     * Performs the genotype mapping per new genotype index.
     *
     * @param newGenotypeIndex the target new genotype index.
     * @param alleleCounts tha correspond to {@code newGenotypeIndex}.
     * @param oldToNewAlleleIndexMap the allele mapping.
     * @param destination where to store the new genotype index mapping to old.
     * @param sortedAlleleCountsBuffer a buffer to re-use to get the genotype-allele-count's sorted allele counts.
     */
    private void genotypeIndexMapPerGenotypeIndex(final int newGenotypeIndex,
                                                  final GenotypeAlleleCounts alleleCounts,
                                                  final int[] oldToNewAlleleIndexMap,
                                                  final int[] destination,
                                                  final int[] sortedAlleleCountsBuffer) {
        final int distinctAlleleCount = alleleCounts.distinctAlleleCount();
        alleleCounts.copyAlleleCounts(sortedAlleleCountsBuffer,0);
        for (int j = 0, jj = 0; j < distinctAlleleCount; j++) {
            final int oldIndex = sortedAlleleCountsBuffer[jj++];
            final int repeats = sortedAlleleCountsBuffer[jj++];
            final int newIndex = oldToNewAlleleIndexMap[oldIndex];
            if (newIndex < 0 || newIndex >= alleleCount)
                throw new IllegalArgumentException("found invalid new allele index (" + newIndex + ") for old index (" + oldIndex + ")");
            for (int k = 0; k < repeats; k++)
                alleleHeap.add(newIndex);
        }
        final int genotypeIndex = alleleHeapToIndex(); // this cleans the heap for the next use.
        destination[newGenotypeIndex] = genotypeIndex;
    }

}
