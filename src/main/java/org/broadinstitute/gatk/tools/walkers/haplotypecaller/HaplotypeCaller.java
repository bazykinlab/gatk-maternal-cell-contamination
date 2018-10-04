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

package org.broadinstitute.gatk.tools.walkers.haplotypecaller;

import com.google.java.contract.Ensures;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.utils.GeneralUtils;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.arguments.DbsnpArgumentCollection;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.engine.filters.BadMateFilter;
import org.broadinstitute.gatk.engine.io.DirectOutputTracker;
import org.broadinstitute.gatk.engine.io.stubs.SAMFileWriterStub;
import org.broadinstitute.gatk.engine.io.stubs.VariantContextWriterStub;
import org.broadinstitute.gatk.engine.iterators.ReadTransformer;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.StandardHCAnnotation;
import org.broadinstitute.gatk.tools.walkers.genotyper.*;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.FixedAFCalculatorProvider;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.gatk.tools.walkers.phasing.PhaseByTransmission;
import org.broadinstitute.gatk.tools.walkers.variantutils.PosteriorLikelihoodsUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.activeregion.ActiveRegion;
import org.broadinstitute.gatk.utils.activeregion.ActiveRegionReadState;
import org.broadinstitute.gatk.utils.activeregion.ActivityProfileState;
import org.broadinstitute.gatk.utils.clipping.ReadClipper;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.AlignmentContextUtils;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.downsampling.DownsamplingUtils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.fragments.FragmentCollection;
import org.broadinstitute.gatk.utils.fragments.FragmentUtils;
import org.broadinstitute.gatk.utils.genotyper.*;
import org.broadinstitute.gatk.utils.gga.GenotypingGivenAllelesUtils;
import org.broadinstitute.gatk.utils.gvcf.GVCFWriter;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.haplotypeBAMWriter.DroppedReadsTracker;
import org.broadinstitute.gatk.utils.haplotypeBAMWriter.HaplotypeBAMWriter;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pairhmm.PairHMM;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.sam.AlignmentUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.gatk.utils.variant.HomoSapiensConstants;

import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

import static org.broadinstitute.gatk.utils.variant.GATKVCFConstants.GENOTYPE_PRIOR_KEY;
import static org.broadinstitute.gatk.utils.variant.GATKVCFConstants.JOINT_LIKELIHOOD_TAG_NAME;
import static org.broadinstitute.gatk.utils.variant.GATKVCFConstants.TRANSMISSION_PROBABILITY_KEY;

/**
 * Call germline SNPs and indels via local re-assembly of haplotypes
 *
 * <p>The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other. It also makes the HaplotypeCaller much better at calling indels than position-based callers like UnifiedGenotyper.</p>

<p>In the so-called GVCF mode used for scalable variant calling in DNA sequence data, HaplotypeCaller runs per-sample to generate an intermediate genomic gVCF (gVCF), which can then be used for joint genotyping of multiple samples in a very efficient way, which enables rapid incremental processing of samples as they roll off the sequencer, as well as scaling to very large cohort sizes (e.g. the 92K exomes of ExAC).</p>

 <p>In addition, HaplotypeCaller is able to handle non-diploid organisms as well as pooled experiment data. Note however that the algorithms used to calculate variant likelihoods is not well suited to extreme allele frequencies (relative to ploidy) so its use is not recommended for somatic (cancer) variant discovery. For that purpose, use MuTect2 instead.</p>

 <p>Finally, HaplotypeCaller is also able to correctly handle the splice junctions that make RNAseq a challenge for most variant callers.</p>
 *
 * <h3>How HaplotypeCaller works</h3>
 *
 * <br />
 * <h4>1. Define active regions </h4>
 *
 * <p>The program determines which regions of the genome it needs to operate on, based on the presence of significant
 * evidence for variation.</p>
 *
 * <br />
 * <h4>2. Determine haplotypes by assembly of the active region </h4>
 *
 * <p>For each ActiveRegion, the program builds a De Bruijn-like graph to reassemble the ActiveRegion, and identifies
 * what are the possible haplotypes present in the data. The program then realigns each haplotype against the reference
 * haplotype using the Smith-Waterman algorithm in order to identify potentially variant sites. </p>
 *
 * <br />
 * <h4>3. Determine likelihoods of the haplotypes given the read data </h4>
 *
 * <p>For each ActiveRegion, the program performs a pairwise alignment of each read against each haplotype using the
 * PairHMM algorithm. This produces a matrix of likelihoods of haplotypes given the read data. These likelihoods are
 * then marginalized to obtain the likelihoods of alleles for each potentially variant site given the read data.   </p>
 *
 * <br />
 * <h4>4. Assign sample genotypes </h4>
 *
 * <p>For each potentially variant site, the program applies Bayes' rule, using the likelihoods of alleles given the
 * read data to calculate the likelihoods of each genotype per sample given the read data observed for that
 * sample. The most likely genotype is then assigned to the sample.    </p>
 *
 * <h3>Input</h3>
 * <p>
 * Input bam file(s) from which to make calls
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * Either a VCF or gVCF file with raw, unfiltered SNP and indel calls. Regular VCFs must be filtered either by variant
 * recalibration (best) or hard-filtering before use in downstream analyses. If using the reference-confidence model
 * workflow for cohort analysis, the output is a GVCF file that must first be run through GenotypeGVCFs and then
 * filtering before further analysis.
 * </p>
 *
 * <h3>Usage examples</h3>
 *
 * <p>These are example commands that show how to run HaplotypeCaller for typical use cases. Square brackets ("[ ]")
 * indicate optional arguments. Note that parameter values shown here may not be the latest recommended; see the
 * Best Practices documentation for detailed recommendations. </p>
 *
 * <br />
 * <h4>Single-sample GVCF calling on DNAseq (for `-ERC GVCF` cohort analysis workflow)</h4>
 * <pre>
 *   java -jar GenomeAnalysisTK.jar \
 *     -R reference.fasta \
 *     -T HaplotypeCaller \
 *     -I sample1.bam \
 *     --emitRefConfidence GVCF \
 *     [--dbsnp dbSNP.vcf] \
 *     [-L targets.interval_list] \
 *     -o output.raw.snps.indels.g.vcf
 * </pre>
 *
 * <h4>Single-sample GVCF calling on DNAseq with allele-specific annotations (for allele-specific cohort analysis workflow)</h4>
 * <pre>
 *   java -jar GenomeAnalysisTK.jar \
 *     -R reference.fasta \
 *     -T HaplotypeCaller \
 *     -I sample1.bam \
 *     --emitRefConfidence GVCF \
 *     [--dbsnp dbSNP.vcf] \
 *     [-L targets.interval_list] \
 *     -G Standard -G AS_Standard \
 *     -o output.raw.snps.indels.AS.g.vcf
 * </pre>
 *
 * <h4>Variant-only calling on DNAseq</h4>
 * <pre>
 *   java -jar GenomeAnalysisTK.jar \
 *     -R reference.fasta \
 *     -T HaplotypeCaller \
 *     -I sample1.bam [-I sample2.bam ...] \
 *     [--dbsnp dbSNP.vcf] \
 *     [-stand_call_conf 30] \
 *     [-L targets.interval_list] \
 *     -o output.raw.snps.indels.vcf
 * </pre>
 *
 * <h4>Variant-only calling on RNAseq</h4>
 * <pre>
 *   java -jar GenomeAnalysisTK.jar \
 *     -R reference.fasta \
 *     -T HaplotypeCaller \
 *     -I sample1.bam \
 *     [--dbsnp dbSNP.vcf] \
 *     -stand_call_conf 20 \
 *     -o output.raw.snps.indels.vcf
 * </pre>
 *
 * <h3>Caveats</h3>
 * <ul>
 * <li>We have not yet fully tested the interaction between the GVCF-based calling or the multisample calling and the
 * RNAseq-specific functionalities. Use those in combination at your own risk.</li>
 * <li>Many users have reported issues running HaplotypeCaller with the -nct argument, so we recommend using Queue to
 * parallelize HaplotypeCaller instead of multithreading.</li>
 * </ul>
 *
 * <h3>Special note on ploidy</h3>
 * <p>This tool is able to handle almost any ploidy (except very high ploidies in large pooled experiments); the ploidy can be specified using the -ploidy argument for non-diploid organisms.</p>
 *
 * <h3>Additional Notes</h3>
 * <ul>
 *     <li>When working with PCR-free data, be sure to set `-pcr_indel_model NONE` (see argument below).</li>
 *     <li>When running in `-ERC GVCF` or `-ERC BP_RESOLUTION` modes, the emitting and calling confidence thresholds
 *     are automatically set to 0. This cannot be overridden by the command line. The thresholds can be set manually
 *     to the desired levels in the next step of the workflow (GenotypeGVCFs)</li>
 * </ul>
 *
 */

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARDISC, extraDocs = {CommandLineGATK.class} )
@PartitionBy(PartitionType.LOCUS)
@BAQMode(ApplicationTime = ReadTransformer.ApplicationTime.FORBIDDEN)
@ActiveRegionTraversalParameters(extension=100, maxRegion=300)
@ReadFilters({HCMappingQualityFilter.class})
@Downsample(by= DownsampleType.BY_SAMPLE, toCoverage=500)
public class HaplotypeCaller extends ActiveRegionWalker<List<VariantContext>, Integer> implements AnnotatorCompatible, NanoSchedulable {
    // -----------------------------------------------------------------------------------------------
    // general haplotype caller arguments
    // -----------------------------------------------------------------------------------------------

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    /**
     * A raw, unfiltered, highly sensitive callset in VCF format.
     */
    @Output(doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;

    @Hidden
    @Advanced
    @Argument(fullName="likelihoodCalculationEngine",shortName="likelihoodEngine",
            doc="What likelihood calculation engine to use to calculate the relative likelihood of reads vs haplotypes",required=false)
    protected ReadLikelihoodCalculationEngine.Implementation likelihoodEngineImplementation = ReadLikelihoodCalculationEngine.Implementation.PairHMM;

    @Hidden
    @Advanced
    @Argument(fullName="heterogeneousKmerSizeResolution",shortName="hksr",doc="How to solve heterogeneous kmer situations using the fast method",required=false)
    protected HeterogeneousKmerSizeResolution heterogeneousKmerSizeResolution = HeterogeneousKmerSizeResolution.COMBO_MIN;

    private HaplotypeBAMWriter haplotypeBAMWriter;

    /**
     * rsIDs from this file are used to populate the ID column of the output. Also, the DB INFO flag will be set when appropriate.
     * dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();
    private double log10GlobalReadMismappingRate;

    /**
     * Active region trimmer reference.
     */
    @ArgumentCollection
    protected ActiveRegionTrimmer trimmer = new ActiveRegionTrimmer();

    public RodBinding<VariantContext> getDbsnpRodBinding() { return dbsnp.dbsnp; }

    /**
     * If a call overlaps with a record from the provided comp track, the INFO field will be annotated
     * as such in the output with the track name (e.g. -comp:FOO will have 'FOO' in the INFO field). Records that are
     * filtered in the comp track will be ignored. Note that 'dbSNP' has been special-cased (see the --dbsnp argument).
     */
    @Advanced
    @Input(fullName="comp", shortName = "comp", doc="Comparison VCF file", required=false)
    public List<RodBinding<VariantContext>> comps = Collections.emptyList();
    public List<RodBinding<VariantContext>> getCompRodBindings() { return comps; }

    // The following are not used by the Unified Genotyper
    public RodBinding<VariantContext> getSnpEffRodBinding() { return null; }
    public List<RodBinding<VariantContext>> getResourceRodBindings() { return Collections.emptyList(); }
    public boolean alwaysAppendDbsnpId() { return false; }

    /**
     * Which annotations to add to the output VCF file. The single value 'none' removes the default annotations.
     * See the VariantAnnotator -list argument to view available annotations.
     */
    @Advanced
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", required=false)
    protected List<String> annotationsToUse = new ArrayList<>();

    /**
     * Which annotations to exclude from output in the VCF file.  Note that this argument has higher priority than the
     * -A or -G arguments, so these annotations will be excluded even if they are explicitly included with the other
     * options. When HaplotypeCaller is run with -ERC GVCF or -ERC BP_RESOLUTION, some annotations are excluded from the
     * output by default because they will only be meaningful once they have been recalculated by GenotypeGVCFs. As
     * of version 3.3 this concerns ChromosomeCounts, FisherStrand, StrandOddsRatio and QualByDepth.
     *
     */
    @Advanced
    @Argument(fullName="excludeAnnotation", shortName="XA", doc="One or more specific annotations to exclude", required=false)
    protected List<String> annotationsToExclude = new ArrayList<>();

    /**
     * Which groups of annotations to add to the output VCF file. The single value 'none' removes the default group. See
     * the VariantAnnotator -list argument to view available groups. Note that this usage is not recommended because
     * it obscures the specific requirements of individual annotations. Any requirements that are not met (e.g. failing
     * to provide a pedigree file for a pedigree-based annotation) may cause the run to fail.
     */
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
    protected List<String> annotationGroupsToUse = new ArrayList<>(Arrays.asList(new String[]{StandardAnnotation.class.getSimpleName(), StandardHCAnnotation.class.getSimpleName() }));

    @ArgumentCollection
    private HaplotypeCallerArgumentCollection HCAC = new HaplotypeCallerArgumentCollection();

    @ArgumentCollection
    private LikelihoodEngineArgumentCollection LEAC = new LikelihoodEngineArgumentCollection();

    /**
     * You can use this argument to specify that HC should process a single sample out of a multisample BAM file. This
     * is especially useful if your samples are all in the same file but you need to run them individually through HC
     * in -ERC GVC mode (which is the recommended usage). Note that the name is case-sensitive.
     */
    @Argument(fullName="sample_name", shortName = "sn", doc="Name of single sample to use from a multi-sample bam", required=false)
    protected String sampleNameToUse = null;

    @ArgumentCollection
    private ReadThreadingAssemblerArgumentCollection RTAC = new ReadThreadingAssemblerArgumentCollection();

    // -----------------------------------------------------------------------------------------------
    // general advanced arguments to control haplotype caller behavior
    // -----------------------------------------------------------------------------------------------


    /**
     * When HC is run in reference confidence mode with banding compression enabled (-ERC GVCF), homozygous-reference
     * sites are compressed into bands of similar genotype quality (GQ) that are emitted as a single VCF record. See
     * the FAQ documentation for more details about the GVCF format.
     *
     * This argument allows you to set the GQ bands. HC expects a list of strictly increasing GQ values
     * that will act as exclusive upper bounds for the GQ bands. To pass multiple values,
     * you provide them one by one with the argument, as in `-GQB 10 -GQB 20 -GQB 30` and so on
     * (this would set the GQ bands to be `[0, 10), [10, 20), [20, 30)` and so on, for example).
     * Note that GQ values are capped at 99 in the GATK, so values must be integers in [1, 100].
     * If the last value is strictly less than 100, the last GQ band will start at that value (inclusive)
     * and end at 100 (exclusive).
     */
    @Advanced
    @Argument(fullName="GVCFGQBands", shortName="GQB", doc="Exclusive upper bounds for reference confidence GQ bands " +
            "(must be in [1, 100] and specified in increasing order)", required = false)
    protected List<Integer> GVCFGQBands = new ArrayList<Integer>(70) {{
        for (int i=1; i<=60; ++i) add(i);
        add(70); add(80); add(90); add(99);
    }};

    /**
     * This parameter determines the maximum size of an indel considered as potentially segregating in the
     * reference model.  It is used to eliminate reads from being indel informative at a site, and determines
     * by that mechanism the certainty in the reference base.  Conceptually, setting this parameter to
     * X means that each informative read is consistent with any indel of size < X being present at a specific
     * position in the genome, given its alignment to the reference.
     */
    @Advanced
    @Argument(fullName="indelSizeToEliminateInRefModel", shortName="ERCIS", doc="The size of an indel to check for in the reference model", required = false)
    protected int indelSizeToEliminateInRefModel = 10;

    // -----------------------------------------------------------------------------------------------
    // general advanced arguments to control haplotype caller behavior
    // -----------------------------------------------------------------------------------------------

    /**
     * Bases with a quality below this threshold will not be used for calling.
     */
    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for calling", required = false)
    public byte MIN_BASE_QUALTY_SCORE = 10;


    /**
     * If this flag is provided, the HaplotypeCaller will include unmapped reads (that have chromosomal coordinates) in the assembly and calling
     * when these reads occur in the region being analyzed.  This situation can occur in paired end analyses, when one read in the read pair
     * gets mapped but its mate is too divergent. In that case, the mate will be marked as unmapped and placed next to the first read, assigned to the same
     * contig and alignment start.  If this flag is provided, the HaplotypeCaller will see such reads, and may make use of them in assembly and calling, where possible.
     */
    @Hidden
    @Argument(fullName="includeUmappedReads", shortName="unmapped", doc="Include unmapped reads with chromosomal coordinates", required = false)
    protected boolean includeUnmappedReads = false;

    @Advanced
    @Argument(fullName="useAllelesTrigger", shortName="allelesTrigger", doc = "Use additional trigger on variants found in an external alleles file", required=false)
    protected boolean USE_ALLELES_TRIGGER = false;

    /**
     * As of GATK 3.3, HaplotypeCaller outputs physical (read-based) information (see version 3.3 release notes and documentation for details). This argument disables that behavior.
     */
    @Advanced
    @Argument(fullName="doNotRunPhysicalPhasing", shortName="doNotRunPhysicalPhasing", doc="Disable physical phasing", required = false)
    protected boolean doNotRunPhysicalPhasing = false;

    // -----------------------------------------------------------------------------------------------
    // arguments for debugging / developing the haplotype caller
    // -----------------------------------------------------------------------------------------------

    @Hidden
    @Argument(fullName="keepRG", shortName="keepRG", doc="Only use reads from this read group when making calls (but use all reads to build the assembly)", required = false)
    protected String keepRG = null;

    /**
     * This argument is intended for benchmarking and scalability testing.
     */
    @Hidden
    @Argument(fullName="justDetermineActiveRegions", shortName="justDetermineActiveRegions", doc = "Just determine ActiveRegions, don't perform assembly or calling", required=false)
    protected boolean justDetermineActiveRegions = false;

    /**
     * This argument is intended for benchmarking and scalability testing.
     */
    @Hidden
    @Argument(fullName="dontGenotype", shortName="dontGenotype", doc = "Perform assembly but do not genotype variants", required=false)
    protected boolean dontGenotype = false;

    @Advanced
    @Argument(fullName="dontUseSoftClippedBases", shortName="dontUseSoftClippedBases", doc="Do not analyze soft clipped bases in the reads", required = false)
    protected boolean dontUseSoftClippedBases = false;

    @Hidden
    @Argument(fullName="captureAssemblyFailureBAM", shortName="captureAssemblyFailureBAM", doc="Write a BAM called assemblyFailure.bam capturing all of the reads that were in the active region when the assembler failed for any reason", required = false)
    protected boolean captureAssemblyFailureBAM = false;

    // Parameters to control read error correction
    /**
     * Enabling this argument may cause fundamental problems with the assembly graph itself.
     */
    @Hidden
    @Argument(fullName="errorCorrectReads", shortName="errorCorrectReads", doc = "Use an exploratory algorithm to error correct the kmers used during assembly", required=false)
    protected boolean errorCorrectReads = false;

    /**
     * When calculating the likelihood of variants, we can try to correct for PCR errors that cause indel artifacts.
     * The correction is based on the reference context, and acts specifically around repetitive sequences that tend
     * to cause PCR errors). The variant likelihoods are penalized in increasing scale as the context around a
     * putative indel is more repetitive (e.g. long homopolymer). The correction can be disabling by specifying
     * '-pcrModel NONE'; in that case the default base insertion/deletion qualities will be used (or taken from the
     * read if generated through the BaseRecalibrator). <b>VERY IMPORTANT: when using PCR-free sequencing data we
     * definitely recommend setting this argument to NONE</b>.
     */
    @Advanced
    @Argument(fullName = "pcr_indel_model", shortName = "pcrModel", doc = "The PCR indel model to use", required = false)
    public PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL pcrErrorModel = PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL.CONSERVATIVE;

    // -----------------------------------------------------------------------------------------------
    // done with Haplotype caller parameters
    // -----------------------------------------------------------------------------------------------

    // the UG engines
    private UnifiedGenotypingEngine activeRegionEvaluationGenotyperEngine = null;

    // the assembly engine
    private LocalAssemblyEngine assemblyEngine = null;

    // the likelihoods engine
    private ReadLikelihoodCalculationEngine likelihoodCalculationEngine = null;

    // the genotyping engine
    private HaplotypeCallerGenotypingEngine genotypingEngine = null;

    // fasta reference reader to supplement the edges of the reference sequence
    protected ReferenceSequenceFile referenceReader;

    // reference base padding size
    private static final int REFERENCE_PADDING = 500;

    /**
     * When downsampling, level the coverage of the reads in each sample to no more than maxReadsInRegionPerSample reads,
     * not reducing coverage at any read start to less than minReadsPerAlignmentStart
     */
    @Argument(fullName = "maxReadsInRegionPerSample", shortName = "maxReadsInRegionPerSample", doc="Maximum reads in an active region", required = false)
    protected int maxReadsInRegionPerSample = 10000;

    @Argument(fullName = "minReadsPerAlignmentStart", shortName = "minReadsPerAlignStart", doc="Minimum number of reads sharing the same alignment start for each genomic location in an active region", required = false)
    protected int minReadsPerAlignmentStart = 10;

    private byte MIN_TAIL_QUALITY;
    private static final byte MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION = 6;

    /**
     * Minimum (exclusive) average number of high quality bases per soft-clip to consider that a set of soft-clips is a
     * high quality set.
     */
    private static final double AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD = 6.0;

    /**
     * Maximum-mininum confidence on a variant to exist to consider the position as a potential variant harbouring locus
     * when looking for active regions.
     */
    private static final double MAXMIN_CONFIDENCE_FOR_CONSIDERING_A_SITE_AS_POSSIBLE_VARIANT_IN_ACTIVE_REGION_DISCOVERY = 4.0;

    /**
     * Minimum ploidy assumed when looking for loci that may harbour variation to identify active regions.
     * <p>
     * By default we take the putative ploidy provided by the user, but this turned out to be too insensitive
     * for low ploidy, notoriously with haploid samples. Therefore we impose this minimum.
     * </p>
     */
    private static final int MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY = 2;

    /**
     * Reads with length lower than this number, after clipping off overhands outside the active region,
     * won't be considered for genotyping.
     */
    private final static int READ_LENGTH_FILTER_THRESHOLD = 10;

    /**
     * Reads with mapping quality lower than this number won't be considered for genotyping, even if they have
     * passed earlier filters (e.g. the walkers' own min MQ filter).
     */
    private static final int READ_QUALITY_FILTER_THRESHOLD = 20;

    private SampleList samplesList;

    private final static Allele FAKE_REF_ALLELE = Allele.create("N", true); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file
    private final static Allele FAKE_ALT_ALLELE = Allele.create("<FAKE_ALT>", false); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file

    ReferenceConfidenceModel referenceConfidenceModel = null;

    /*
     * Is the bad mate filter disabled by the argument -drf BadMate?
     */
    private boolean isBadMateFilterDisabled = false;


    //Matrix of priors for all genotype combinations
    private EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>> mvCountMatrix;

    //Matrix of allele transmission
    private EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType,TrioPhase>>> transmissionMatrix;

    ////////////////////////////////////////////////////////////////////////
    //// Deprecated Arguments					                        ////
    //// Keeping them here is meant to provide informative error messages //
    //// when an argument has been put out of service		            ////
    ////////////////////////////////////////////////////////////////////////
    /**
     * @deprecated
     * Deprecated: 2015-04-01, J.White
     * mergeVariantsViaLD = false made final
     */
    @Hidden
    @Deprecated
    @Argument(fullName="mergeVariantsViaLD", shortName="mergeVariantsViaLD", doc="DEPRECATED; This argument is no longer used in GATK versions 3.4 and newer. Please see the online documentation for the latest usage recommendations.", required = false)
    static final boolean mergeVariantsViaLD = false;

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------
    public final double NO_TRANSMISSION_PROB = -1.0;
    private enum FamilyMember {
        MOTHER,
        FATHER,
        CHILD
    }
    private double deNovoPrior=1e-8;
    private boolean fatherFAlleleFirst=false;
    //Random number generator
    private Random rand = new Random();

    private int getCombinationMVCount(GenotypeType mother, GenotypeType father, GenotypeType child){

        //Child is no call => No MV
        if(child == GenotypeType.NO_CALL || child == GenotypeType.UNAVAILABLE)
            return 0;
        //Add parents with genotypes for the evaluation
        ArrayList<GenotypeType> parents = new ArrayList<GenotypeType>();
        if (!(mother == GenotypeType.NO_CALL || mother == GenotypeType.UNAVAILABLE))
            parents.add(mother);
        if (!(father == GenotypeType.NO_CALL || father == GenotypeType.UNAVAILABLE))
            parents.add(father);

        //Both parents no calls => No MV
        if (parents.isEmpty())
            return 0;

        //If at least one parent had a genotype, then count the number of ref and alt alleles that can be passed
        int parentsNumRefAlleles = 0;
        int parentsNumAltAlleles = 0;

        for(GenotypeType parent : parents){
            if(parent == GenotypeType.HOM_REF){
                parentsNumRefAlleles++;
            }
            else if(parent == GenotypeType.HET){
                parentsNumRefAlleles++;
                parentsNumAltAlleles++;
            }
            else if(parent == GenotypeType.HOM_VAR){
                parentsNumAltAlleles++;
            }
        }

        //Case Child is HomRef
        if(child == GenotypeType.HOM_REF){
            if(parentsNumRefAlleles == parents.size())
                return 0;
            else return (parents.size()-parentsNumRefAlleles);
        }

        //Case child is HomVar
        if(child == GenotypeType.HOM_VAR){
            if(parentsNumAltAlleles == parents.size())
                return 0;
            else return parents.size()-parentsNumAltAlleles;
        }

        //Case child is Het
        if(child == GenotypeType.HET && ((parentsNumRefAlleles > 0 && parentsNumAltAlleles > 0) || parents.size()<2))
            return 0;

        //MV
        return 1;
    }
    //Stores a conceptual trio or parent/child pair genotype combination along with its phasing.
    //This combination can then be "applied" to a given trio or pair using the getPhasedGenotypes method.
    public class TrioPhase {

        //Create 2 fake alleles
        //The actual bases will never be used but the Genotypes created using the alleles will be.
        private final Allele REF = Allele.create("A",true);
        private final Allele VAR = Allele.create("A",false);
        private final Allele NO_CALL = Allele.create(".",false);
        private final String DUMMY_NAME = "DummySample";

        private EnumMap<FamilyMember,Genotype> trioPhasedGenotypes = new EnumMap<FamilyMember, Genotype>(FamilyMember.class);

        private ArrayList<Allele> getAlleles(GenotypeType genotype){
            ArrayList<Allele> alleles = new ArrayList<Allele>(2);
            if(genotype == GenotypeType.HOM_REF){
                alleles.add(REF);
                alleles.add(REF);
            }
            else if(genotype == GenotypeType.HET){
                alleles.add(REF);
                alleles.add(VAR);
            }
            else if(genotype == GenotypeType.HOM_VAR){
                alleles.add(VAR);
                alleles.add(VAR);
            }
            else{
                return null;
            }
            return alleles;
        }

        private boolean isPhasable(GenotypeType genotype){
            return genotype == GenotypeType.HOM_REF || genotype == GenotypeType.HET || genotype == GenotypeType.HOM_VAR;
        }

        //Create a new Genotype based on information from a single individual
        //Homozygous genotypes will be set as phased, heterozygous won't be
        private void phaseSingleIndividualAlleles(GenotypeType genotype, FamilyMember familyMember){
            boolean phase = genotype == GenotypeType.HOM_REF || genotype == GenotypeType.HOM_VAR;
            trioPhasedGenotypes.put(familyMember, makeGenotype(genotype, phase));
        }

        private Genotype makeGenotype(final GenotypeType type, boolean phase) {
            return makeGenotype(getAlleles(type), phase);
        }

        private Genotype makeGenotype(final List<Allele> alleles, boolean phase) {
            final GenotypeBuilder gb = new GenotypeBuilder(DUMMY_NAME, alleles);
            gb.phased(phase);
            return gb.make();
        }

        //Find the phase for a parent/child pair
        private void phasePairAlleles(GenotypeType parentGenotype, GenotypeType childGenotype, FamilyMember parent){

            //Special case for Het/Het as it is ambiguous
            if(parentGenotype == GenotypeType.HET && childGenotype == GenotypeType.HET){
                trioPhasedGenotypes.put(parent, makeGenotype(parentGenotype, false));
                trioPhasedGenotypes.put(FamilyMember.CHILD, makeGenotype(childGenotype, false));
                return;
            }

            ArrayList<Allele> parentAlleles = getAlleles(parentGenotype);
            ArrayList<Allele> childAlleles = getAlleles(childGenotype);
            ArrayList<Allele> parentPhasedAlleles = new ArrayList<Allele>(2);
            ArrayList<Allele> childPhasedAlleles = new ArrayList<Allele>(2);

            //If there is a possible phasing between the parent and child => phase
            int childTransmittedAlleleIndex = childAlleles.indexOf(parentAlleles.get(0));
            if(childTransmittedAlleleIndex > -1){
                trioPhasedGenotypes.put(parent, makeGenotype(parentAlleles, true));
                childPhasedAlleles.add(childAlleles.remove(childTransmittedAlleleIndex));
                if(parent.equals(FamilyMember.MOTHER))
                    childPhasedAlleles.add(childAlleles.get(0));
                else
                    childPhasedAlleles.add(0,childAlleles.get(0));
                trioPhasedGenotypes.put(FamilyMember.CHILD, makeGenotype(childPhasedAlleles, true));
            }
            else if((childTransmittedAlleleIndex = childAlleles.indexOf(parentAlleles.get(1))) > -1){
                parentPhasedAlleles.add(parentAlleles.get(1));
                parentPhasedAlleles.add(parentAlleles.get(0));
                trioPhasedGenotypes.put(parent, makeGenotype(parentPhasedAlleles, true));
                childPhasedAlleles.add(childAlleles.remove(childTransmittedAlleleIndex));
                if(parent.equals(FamilyMember.MOTHER))
                    childPhasedAlleles.add(childAlleles.get(0));
                else
                    childPhasedAlleles.add(0,childAlleles.get(0));
                trioPhasedGenotypes.put(FamilyMember.CHILD, makeGenotype(childPhasedAlleles, true));
            }
            //This is a Mendelian Violation => Do not phase
            else{
                trioPhasedGenotypes.put(parent, makeGenotype(parentGenotype, false));
                trioPhasedGenotypes.put(FamilyMember.CHILD, makeGenotype(childGenotype, false));
            }
        }

        //Phases a family by transmission
        private void phaseFamilyAlleles(GenotypeType mother, GenotypeType father, GenotypeType child){

            Set<ArrayList<Allele>> possiblePhasedChildGenotypes = new HashSet<ArrayList<Allele>>();
            ArrayList<Allele> motherAlleles = getAlleles(mother);
            ArrayList<Allele> fatherAlleles = getAlleles(father);
            ArrayList<Allele> childAlleles = getAlleles(child);

            //Build all possible child genotypes for the given parent's genotypes
            for (Allele momAllele : motherAlleles) {
                for (Allele fatherAllele : fatherAlleles) {
                    ArrayList<Allele> possiblePhasedChildAlleles = new ArrayList<Allele>(2);
                    possiblePhasedChildAlleles.add(momAllele);
                    possiblePhasedChildAlleles.add(fatherAllele);
                    possiblePhasedChildGenotypes.add(possiblePhasedChildAlleles);
                }
            }

            for (ArrayList<Allele> childPhasedAllelesAlleles : possiblePhasedChildGenotypes) {
                int firstAlleleIndex = childPhasedAllelesAlleles.indexOf(childAlleles.get(0));
                int secondAlleleIndex = childPhasedAllelesAlleles.lastIndexOf(childAlleles.get(1));
                //If a possible combination has been found, create the genotypes
                if (firstAlleleIndex != secondAlleleIndex && firstAlleleIndex > -1 && secondAlleleIndex > -1) {
                    //Create mother's genotype
                    ArrayList<Allele> motherPhasedAlleles = new ArrayList<Allele>(2);
                    motherPhasedAlleles.add(childPhasedAllelesAlleles.get(0));
                    if(motherAlleles.get(0) != motherPhasedAlleles.get(0))
                        motherPhasedAlleles.add(motherAlleles.get(0));
                    else
                        motherPhasedAlleles.add(motherAlleles.get(1));
                    trioPhasedGenotypes.put(FamilyMember.MOTHER, makeGenotype(motherPhasedAlleles, true));

                    //Create father's genotype
                    ArrayList<Allele> fatherPhasedAlleles = new ArrayList<Allele>(2);
                    fatherPhasedAlleles.add(childPhasedAllelesAlleles.get(1));
                    if(fatherAlleles.get(0) != fatherPhasedAlleles.get(0))
                        fatherPhasedAlleles.add(fatherAlleles.get(0));
                    else
                        fatherPhasedAlleles.add(fatherAlleles.get(1));
                    trioPhasedGenotypes.put(FamilyMember.FATHER, makeGenotype(fatherPhasedAlleles,true));

                    //Create child's genotype
                    trioPhasedGenotypes.put(FamilyMember.CHILD, makeGenotype(childPhasedAllelesAlleles,true));

                    //Once a phased combination is found; exit
                    return;
                }
            }

            //If this is reached then no phasing could be found
            trioPhasedGenotypes.put(FamilyMember.MOTHER, makeGenotype(mother,false));
            trioPhasedGenotypes.put(FamilyMember.FATHER, makeGenotype(father,false));
            trioPhasedGenotypes.put(FamilyMember.CHILD, makeGenotype(child,false));
        }

        /*  Constructor: Creates a conceptual trio genotype combination from the given genotypes.
            If one or more genotypes are set as NO_CALL or UNAVAILABLE, it will phase them like a pair
            or single individual.
        */
        public TrioPhase(GenotypeType mother, GenotypeType father, GenotypeType child){

            //Take care of cases where one or more family members are no call
            if(!isPhasable(child)){

                phaseSingleIndividualAlleles(mother, FamilyMember.MOTHER);
                phaseSingleIndividualAlleles(father, FamilyMember.FATHER);
                phaseSingleIndividualAlleles(child, FamilyMember.CHILD);
            }
            else if(!isPhasable(mother)){

                phaseSingleIndividualAlleles(mother, FamilyMember.MOTHER);
                if(!isPhasable(father)){
                    phaseSingleIndividualAlleles(father, FamilyMember.FATHER);
                    phaseSingleIndividualAlleles(child, FamilyMember.CHILD);
                }
                else
                    phasePairAlleles(father, child, FamilyMember.FATHER);
            }
            else if(!isPhasable(father)){

                phasePairAlleles(mother, child, FamilyMember.MOTHER);
                phaseSingleIndividualAlleles(father, FamilyMember.FATHER);
            }
            //Special case for Het/Het/Het as it is ambiguous
            else if(mother == GenotypeType.HET && father  == GenotypeType.HET && child == GenotypeType.HET){

                phaseSingleIndividualAlleles(mother, FamilyMember.MOTHER);
                phaseSingleIndividualAlleles(father, FamilyMember.FATHER);
                phaseSingleIndividualAlleles(child, FamilyMember.CHILD);
            }
            //All family members have genotypes and at least one of them is not Het
            else{
                System.out.println("TrioPhase --- phaseFamilyAlleles" );
                phaseFamilyAlleles(mother, father, child);
            }

            //If child should phased genotype should be father first, then swap the alleles
            if(fatherFAlleleFirst && trioPhasedGenotypes.get(FamilyMember.CHILD).isPhased()){
                ArrayList<Allele> childAlleles = new ArrayList<Allele>(trioPhasedGenotypes.get(FamilyMember.CHILD).getAlleles());
                childAlleles.add(childAlleles.remove(0));
                trioPhasedGenotypes.put(FamilyMember.CHILD,makeGenotype(childAlleles,true));
            }

        }

        /**
         * Applies the trio genotype combination to the given trio.
         * @param ref: Reference allele
         * @param alt: Alternate allele
         * @param motherGenotype: Genotype of the mother to phase using this trio genotype combination
         * @param fatherGenotype: Genotype of the father to phase using this trio genotype combination
         * @param childGenotype: Genotype of the child to phase using this trio genotype combination
         * @param transmissionProb: Probability for this trio genotype combination to be correct (pass NO_TRANSMISSION_PROB if unavailable)
         * @param phasedGenotypes: An ArrayList<Genotype> to which the newly phased genotypes are added in the following order: Mother, Father, Child
         */
        public void getPhasedGenotypes(Allele ref, Allele alt, Genotype motherGenotype, Genotype fatherGenotype, Genotype childGenotype, double transmissionProb,ArrayList<Genotype> phasedGenotypes){
            System.out.println("transmissionProb" + transmissionProb);
            phasedGenotypes.add(getPhasedGenotype(ref,alt,motherGenotype,transmissionProb,this.trioPhasedGenotypes.get(FamilyMember.MOTHER)));
            phasedGenotypes.add(getPhasedGenotype(ref,alt,fatherGenotype,transmissionProb,this.trioPhasedGenotypes.get(FamilyMember.FATHER)));
            phasedGenotypes.add(getPhasedGenotype(ref,alt,childGenotype,transmissionProb,this.trioPhasedGenotypes.get(FamilyMember.CHILD)));
        }

        private Genotype getPhasedGenotype(Allele refAllele, Allele altAllele, Genotype genotype, double transmissionProb, Genotype phasedGenotype){
            System.out.println("getPhasedGenotype(" + refAllele + ", " + altAllele + "," + genotype + ", " + transmissionProb + "," + phasedGenotype);

            int phredScoreTransmission = -1;
            if(transmissionProb != NO_TRANSMISSION_PROB){
                double dphredScoreTransmission = Byte.MAX_VALUE;
                try {
                    dphredScoreTransmission = QualityUtils.phredScaleLog10ErrorRate(Math.log10(1 - (transmissionProb)));
                } catch(Exception e ) {
                    System.out.println("EXCEPTION transmissionProb =" + transmissionProb + " " + refAllele + " -- " + altAllele + " " +
                    genotype.getAlleles());
                }
                phredScoreTransmission = dphredScoreTransmission < Byte.MAX_VALUE ? (byte)dphredScoreTransmission : Byte.MAX_VALUE;
            }
            //Handle null, missing and unavailable genotypes
            //Note that only cases where a null/missing/unavailable genotype was passed in the first place can lead to a null/missing/unavailable
            //genotype so it is safe to return the original genotype in this case.
            //In addition, if the phasing confidence is 0, then return the unphased, original genotypes.
            if(phredScoreTransmission ==0 || genotype == null || !isPhasable(genotype.getType()))
                return genotype;

            //Add the transmission probability
            Map<String, Object> genotypeAttributes = new HashMap<String, Object>();
            genotypeAttributes.putAll(genotype.getExtendedAttributes());
            if(transmissionProb>NO_TRANSMISSION_PROB)
                genotypeAttributes.put(GATKVCFConstants.TRANSMISSION_PROBABILITY_KEY, phredScoreTransmission);

            ArrayList<Allele> phasedAlleles = new ArrayList<Allele>(2);
            for(Allele allele : phasedGenotype.getAlleles()){
                if(allele.isReference())
                    phasedAlleles.add(refAllele);
                else if(allele.isNonReference())
                    phasedAlleles.add(altAllele);
                    //At this point there should not be any other alleles left
                else
                    throw new UserException(String.format("BUG: Unexpected allele: %s. Please report.",allele.toString()));

            }

            //Compute the new Log10Error if the genotype is different from the original genotype
            double log10Error;
            if(genotype.getType() == phasedGenotype.getType())
                log10Error = genotype.getLog10PError();
            else
                log10Error = genotype.getLikelihoods().getLog10GQ(phasedGenotype.getType());

            return new GenotypeBuilder(genotype).alleles(phasedAlleles)
                    .log10PError(log10Error)
                    .attributes(genotypeAttributes)
                    .phased(phasedGenotype.isPhased()).make();
        }


    }
    private void buildMatrices(){
        mvCountMatrix = new EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>>(GenotypeType.class);
        transmissionMatrix = new EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType,TrioPhase>>>(GenotypeType.class);
        for(GenotypeType mother : GenotypeType.values()){
            mvCountMatrix.put(mother,new EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>(GenotypeType.class));
            transmissionMatrix.put(mother,new EnumMap<GenotypeType,EnumMap<GenotypeType,TrioPhase>>(GenotypeType.class));
            for(GenotypeType father : GenotypeType.values()){
                mvCountMatrix.get(mother).put(father,new EnumMap<GenotypeType, Integer>(GenotypeType.class));
                transmissionMatrix.get(mother).put(father,new EnumMap<GenotypeType,TrioPhase>(GenotypeType.class));
                for(GenotypeType child : GenotypeType.values()){
                    mvCountMatrix.get(mother).get(father).put(child, getCombinationMVCount(mother, father, child));
                    transmissionMatrix.get(mother).get(father).put(child,new TrioPhase(mother,father,child));
                }
            }
        }
    }


    public void initialize() {
        super.initialize();

        buildMatrices();

        if (HCAC.genotypeArgs.samplePloidy != HomoSapiensConstants.DEFAULT_PLOIDY && !doNotRunPhysicalPhasing) {
            doNotRunPhysicalPhasing = true;
            logger.info("Currently, physical phasing is not available when ploidy is different than " + HomoSapiensConstants.DEFAULT_PLOIDY + "; therefore it won't be performed");
        }

        if (dontGenotype && emitReferenceConfidence())
            throw new UserException("You cannot request gVCF output and 'do not genotype' at the same time");

        if ( emitReferenceConfidence() ) {

            if (HCAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES)
                throw new UserException.BadArgumentValue("ERC/gt_mode","you cannot request reference confidence output and GENOTYPE_GIVEN_ALLELES at the same time");

            HCAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING = -0.0;

            // also, we don't need to output several of the annotations
            annotationsToExclude.add("ChromosomeCounts");
            annotationsToExclude.add("FisherStrand");
            annotationsToExclude.add("StrandOddsRatio");
            annotationsToExclude.add("QualByDepth");

            // but we definitely want certain other ones
            annotationsToUse.add("StrandBiasBySample");
            logger.info("Standard Emitting and Calling confidence set to 0.0 for reference-model confidence output");
            if (!HCAC.annotateAllSitesWithPLs)
                logger.info("All sites annotated with PLs forced to true for reference-model confidence output");
            HCAC.annotateAllSitesWithPLs = true;
        } else if ( ! doNotRunPhysicalPhasing ) {
            doNotRunPhysicalPhasing = true;
            logger.info("Disabling physical phasing, which is supported only for reference-model confidence output");
        }

        final GenomeAnalysisEngine toolkit = getToolkit();
        samplesList = toolkit.getReadSampleList();
        Set<String> sampleSet = SampleListUtils.asSet(samplesList);

        if (sampleNameToUse != null) {
            if (!sampleSet.contains(sampleNameToUse))
                throw new UserException.BadArgumentValue("sample_name", "Specified name does not exist in input bam files");
            if (sampleSet.size() == 1) {
                //No reason to incur performance penalty associated with filtering if they specified the name of the only sample
                sampleNameToUse = null;
            } else {
                samplesList = new IndexedSampleList(sampleNameToUse);
                sampleSet = SampleListUtils.asSet(samplesList);
            }
        }


        // create a UAC but with the exactCallsLog = null, so we only output the log for the HC caller itself, if requested
        final UnifiedArgumentCollection simpleUAC = HCAC.cloneTo(UnifiedArgumentCollection.class);
        simpleUAC.outputMode = OutputMode.EMIT_VARIANTS_ONLY;
        simpleUAC.genotypingOutputMode = GenotypingOutputMode.DISCOVERY;
        simpleUAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING = Math.min(MAXMIN_CONFIDENCE_FOR_CONSIDERING_A_SITE_AS_POSSIBLE_VARIANT_IN_ACTIVE_REGION_DISCOVERY, HCAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING ); // low values used for isActive determination only, default/user-specified values used for actual calling
        simpleUAC.CONTAMINATION_FRACTION = 0.0;
        simpleUAC.CONTAMINATION_FRACTION_FILE = null;
        simpleUAC.exactCallsLog = null;
        // Seems that at least with some test data we can lose genuine haploid variation if we use
        // UGs engine with ploidy == 1
        simpleUAC.genotypeArgs.samplePloidy = Math.max(MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY, HCAC.genotypeArgs.samplePloidy);

        activeRegionEvaluationGenotyperEngine = new UnifiedGenotypingEngine(simpleUAC,
                FixedAFCalculatorProvider.createThreadSafeProvider(getToolkit(),simpleUAC,logger), toolkit);
        activeRegionEvaluationGenotyperEngine.setLogger(logger);

        if( HCAC.CONTAMINATION_FRACTION_FILE != null )
            HCAC.setSampleContamination(AlleleBiasedDownsamplingUtils.loadContaminationFile(HCAC.CONTAMINATION_FRACTION_FILE, HCAC.CONTAMINATION_FRACTION, sampleSet, logger));

        if( HCAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES && RTAC.consensusMode )
            throw new UserException("HaplotypeCaller cannot be run in both GENOTYPE_GIVEN_ALLELES mode and in consensus mode at the same time. Please choose one or the other.");

        final GenomeLocParser genomeLocParser = toolkit.getGenomeLocParser();

        genotypingEngine = new HaplotypeCallerGenotypingEngine(HCAC, samplesList, genomeLocParser, FixedAFCalculatorProvider.createThreadSafeProvider(getToolkit(), HCAC,logger), !doNotRunPhysicalPhasing);
        // initialize the output VCF header
        final VariantAnnotatorEngine annotationEngine = new VariantAnnotatorEngine(annotationGroupsToUse, annotationsToUse, annotationsToExclude, this, getToolkit());

        final Set<VCFHeaderLine> headerInfo = new HashSet<>();

        //initialize the annotations (this is particularly important to turn off RankSumTest dithering in integration tests)
        //do this before we write the header because SnpEff adds to header lines
        annotationEngine.invokeAnnotationInitializationMethods(headerInfo);

        headerInfo.addAll(genotypingEngine.getAppropriateVCFInfoHeaders());
        // all annotation fields from VariantAnnotatorEngine
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());
        // all callers need to add these standard annotation header lines
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.DOWNSAMPLED_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        // all callers need to add these standard FORMAT field header lines
        VCFStandardHeaderLines.addStandardFormatLines(headerInfo, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);

        if ( ! doNotRunPhysicalPhasing ) {
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));
        }

        // FILTER fields are added unconditionally as it's not always 100% certain the circumstances
        // where the filters are used.  For example, in emitting all sites the lowQual field is used
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.LOW_QUAL_FILTER_NAME));

        headerInfo.add(new VCFInfoHeaderLine(GENOTYPE_PRIOR_KEY, VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Genotype Likelihood Prior"));
        headerInfo.add(new VCFFormatHeaderLine(JOINT_LIKELIHOOD_TAG_NAME, 1, VCFHeaderLineType.Integer, "Phred-scaled joint likelihood of the genotype combination (before applying family priors)"));
        headerInfo.addAll(GATKVCFUtils.getHeaderFields(this.getToolkit()));
        headerInfo.add(new VCFFormatHeaderLine(TRANSMISSION_PROBABILITY_KEY, 1, VCFHeaderLineType.Integer, "Phred score of the genotype combination and phase given that the genotypes are correct"));
        headerInfo.add(new VCFHeaderLine("source", "PhaseByTransmission"));


        initializeReferenceConfidenceModel(samplesList, headerInfo);

        vcfWriter.writeHeader(new VCFHeader(headerInfo, sampleSet));

        // fasta reference reader to supplement the edges of the reference sequence
        referenceReader = CachingIndexedFastaSequenceFile.checkAndCreate(getToolkit().getArguments().referenceFile);

        // create and setup the assembler
        assemblyEngine = new ReadThreadingAssembler(RTAC.maxNumHaplotypesInPopulation, RTAC.kmerSizes, RTAC.dontIncreaseKmerSizesForCycles, RTAC.allowNonUniqueKmersInRef, RTAC.numPruningSamples);

        assemblyEngine.setErrorCorrectKmers(RTAC.errorCorrectKmers);
        assemblyEngine.setPruneFactor(RTAC.MIN_PRUNE_FACTOR);
        assemblyEngine.setDebug(HCAC.DEBUG);
        assemblyEngine.setDebugGraphTransformations(RTAC.debugGraphTransformations);
        assemblyEngine.setAllowCyclesInKmerGraphToGeneratePaths(RTAC.allowCyclesInKmerGraphToGeneratePaths);
        assemblyEngine.setRecoverDanglingBranches(!RTAC.doNotRecoverDanglingBranches);
        assemblyEngine.setMinDanglingBranchLength(RTAC.minDanglingBranchLength);
        assemblyEngine.setMinBaseQualityToUseInAssembly(MIN_BASE_QUALTY_SCORE);

        MIN_TAIL_QUALITY = (byte)(MIN_BASE_QUALTY_SCORE - 1);

        if ( RTAC.graphWriter != null ) assemblyEngine.setGraphWriter(RTAC.graphWriter);

        // setup the likelihood calculation engine
        if ( LEAC.phredScaledGlobalReadMismappingRate < 0 ) LEAC.phredScaledGlobalReadMismappingRate = -1;

        // configure the global mismapping rate
        if ( LEAC.phredScaledGlobalReadMismappingRate < 0 ) {
            log10GlobalReadMismappingRate = - Double.MAX_VALUE;
        } else {
            log10GlobalReadMismappingRate = QualityUtils.qualToErrorProbLog10(LEAC.phredScaledGlobalReadMismappingRate);
            logger.info("Using global mismapping rate of " + LEAC.phredScaledGlobalReadMismappingRate + " => " + log10GlobalReadMismappingRate + " in log10 likelihood units");
        }

        //static member function - set number of threads
        PairHMM.setNumberOfThreads(getToolkit().getTotalNumberOfThreads());
        // create our likelihood calculation engine
        likelihoodCalculationEngine = createLikelihoodCalculationEngine();

        final MergeVariantsAcrossHaplotypes variantMerger = new MergeVariantsAcrossHaplotypes();

        genotypingEngine.setCrossHaplotypeEventMerger(variantMerger);

        genotypingEngine.setAnnotationEngine(annotationEngine);

        if ( HCAC.bamWriter != null ) {
            // we currently do not support multi-threaded BAM writing, so exception out
            if ( getToolkit().getTotalNumberOfThreads() > 1 )
                throw new UserException.BadArgumentValue("bamout", "Currently cannot emit a BAM file from the HaplotypeCaller in multi-threaded mode.");
            haplotypeBAMWriter = HaplotypeBAMWriter.create(HCAC.bamWriterType, HCAC.bamWriter, getToolkit().getSAMFileHeader());
        }

        trimmer.initialize(getToolkit().getGenomeLocParser(), HCAC.DEBUG,
                HCAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES,emitReferenceConfidence());

        isBadMateFilterDisabled = toolkit.getArguments().disabledReadFilters.contains("BadMate");
    }

    private void initializeReferenceConfidenceModel(final SampleList samples, final Set<VCFHeaderLine> headerInfo) {
        referenceConfidenceModel = new ReferenceConfidenceModel(getToolkit().getGenomeLocParser(), samples, getToolkit().getSAMFileHeader(), indelSizeToEliminateInRefModel);
        if ( emitReferenceConfidence() ) {
            if ( samples.sampleCount() != 1 )
                throw new UserException.BadArgumentValue("emitRefConfidence", "Can only be used in single sample mode currently. Use the sample_name argument to run on a single sample out of a multi-sample BAM file.");
            headerInfo.addAll(referenceConfidenceModel.getVCFHeaderLines());
            if ( HCAC.emitReferenceConfidence == ReferenceConfidenceMode.GVCF ) {
                // A kluge to enforce the use of this indexing strategy - must set the gVCF indexing values if not a using a gVCF output file .
                // An output gVCF file automatically sets the indexing values because it has the .g.vcf extension.
                if (!GATKVCFUtils.usingGVCFIndexingArguments(getToolkit().getArguments().variant_index_type, getToolkit().getArguments().variant_index_parameter) && !isGVCF()) {
                    throw new UserException.GVCFIndexException(GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
                }

                try {
                    vcfWriter = new GVCFWriter(vcfWriter, GVCFGQBands, HCAC.genotypeArgs.samplePloidy);
                } catch ( final IllegalArgumentException e ) {
                    throw new UserException.BadArgumentValue("GVCFGQBands", e.getMessage());
                }
            }
        }
    }

    /**
     * Instantiates the appropriate likelihood calculation engine.
     *
     * @return never {@code null}.
     */
    private ReadLikelihoodCalculationEngine createLikelihoodCalculationEngine() {
        switch (likelihoodEngineImplementation) {
            case PairHMM:
                return new PairHMMLikelihoodCalculationEngine( (byte) LEAC.gcpHMM, LEAC.pairHMM, LEAC.pairHMMNativeArgs.getPairHMMArgs(), log10GlobalReadMismappingRate, LEAC.noFpga, pcrErrorModel );
            case GraphBased:
                return new GraphBasedLikelihoodCalculationEngine( (byte) LEAC.gcpHMM,log10GlobalReadMismappingRate, heterogeneousKmerSizeResolution, HCAC.DEBUG, RTAC.debugGraphTransformations);
            case Random:
                return new RandomLikelihoodCalculationEngine();
            default:
                //Note: we do not include in the error message list as it is of no grand public interest.
                throw new UserException("Unsupported likelihood calculation engine '" + likelihoodCalculationEngine +
                        "'. Please use one of the following instead: 'PairHMM' or 'GraphBased'.");
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // isActive
    //
    //---------------------------------------------------------------------------------------------------------------

    // enable deletions in the pileup
    @Override
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    // enable non primary and extended reads in the active region
    @Override
    public EnumSet<ActiveRegionReadState> desiredReadStates() {
        if ( includeUnmappedReads )
            throw new UserException.BadArgumentValue("includeUnmappedReads", "is not yet functional");
        else
            return EnumSet.of(
                    ActiveRegionReadState.PRIMARY,
                    ActiveRegionReadState.NONPRIMARY,
                    ActiveRegionReadState.EXTENDED);
    }

    @Override
    @Ensures({"result.isActiveProb >= 0.0", "result.isActiveProb <= 1.0"})
    public ActivityProfileState isActive( final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context ) {

        if( HCAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES ) {
            final VariantContext vcFromAllelesRod = GenotypingGivenAllelesUtils.composeGivenAllelesVariantContextFromRod(tracker, ref.getLocus(), false, logger, HCAC.alleles);
            if( vcFromAllelesRod != null ) {
                return new ActivityProfileState(ref.getLocus(), 1.0);
            }
        }

        if( USE_ALLELES_TRIGGER ) {
            return new ActivityProfileState( ref.getLocus(), tracker.getValues(HCAC.alleles, ref.getLocus()).size() > 0 ? 1.0 : 0.0 );
        }

        if( context == null || context.getBasePileup().isEmpty() )
            // if we don't have any data, just abort early
            return new ActivityProfileState(ref.getLocus(), 0.0);

        final int ploidy = activeRegionEvaluationGenotyperEngine.getConfiguration().genotypeArgs.samplePloidy;
        final List<Allele> noCall = GATKVariantContextUtils.noCallAlleles(ploidy); // used to noCall all genotypes until the exact model is applied
        final Map<String, AlignmentContext> splitContexts = AlignmentContextUtils.splitContextBySampleName(context);
        final GenotypesContext genotypes = GenotypesContext.create(splitContexts.keySet().size());
        final MathUtils.RunningAverage averageHQSoftClips = new MathUtils.RunningAverage();
        final GenotypingModel genotypingModel = genotypingEngine.getGenotypingModel();
        for( final Map.Entry<String, AlignmentContext> sample : splitContexts.entrySet() ) {
            final String sampleName = sample.getKey();
            // The ploidy here is not dictated by the sample but by the simple genotyping-engine used to determine whether regions are active or not.
            final int activeRegionDetectionHackishSamplePloidy = activeRegionEvaluationGenotyperEngine.getConfiguration().genotypeArgs.samplePloidy;
            final double[] genotypeLikelihoods = referenceConfidenceModel.calcGenotypeLikelihoodsOfRefVsAny(activeRegionDetectionHackishSamplePloidy,sample.getValue().getBasePileup(), ref.getBase(), MIN_BASE_QUALTY_SCORE, averageHQSoftClips).genotypeLikelihoods;
            genotypes.add( new GenotypeBuilder(sample.getKey()).alleles(noCall).PL(genotypeLikelihoods).make() );
        }

        final List<Allele> alleles = Arrays.asList(FAKE_REF_ALLELE , FAKE_ALT_ALLELE);
        final double isActiveProb;

        if (genotypes.size() == 1) {
            // Faster implementation avoiding the costly and over complicated Exact AFCalculator machinery:
            // This is the case when doing GVCF output.
            isActiveProb = activeRegionEvaluationGenotyperEngine.calculateSingleSampleRefVsAnyActiveStateProfileValue(genotypes.get(0).getLikelihoods().getAsVector());
        } else {
            final VariantCallContext vcOut = activeRegionEvaluationGenotyperEngine.calculateGenotypes(new VariantContextBuilder("HCisActive!", context.getContig(), context.getLocation().getStart(), context.getLocation().getStop(), alleles).genotypes(genotypes).make(), GenotypeLikelihoodsCalculationModel.Model.SNP);
            isActiveProb = vcOut == null ? 0.0 : QualityUtils.qualToProb(vcOut.getPhredScaledQual());
        }
        return new ActivityProfileState( ref.getLocus(), isActiveProb, averageHQSoftClips.mean() > AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD ? ActivityProfileState.Type.HIGH_QUALITY_SOFT_CLIPS : ActivityProfileState.Type.NONE, averageHQSoftClips.mean() );
    }

    // From PhaseByTransmission

    private EnumMap<GenotypeType,Double> getLikelihoodsAsMapSafeNull(Genotype genotype){
        if(genotype == null || !genotype.isCalled() || genotype.getLikelihoods() == null){
            EnumMap<GenotypeType,Double> likelihoods = new EnumMap<GenotypeType, Double>(GenotypeType.class);
            likelihoods.put(GenotypeType.HOM_REF,1.0/3.0);
            likelihoods.put(GenotypeType.HET,1.0/3.0);
            likelihoods.put(GenotypeType.HOM_VAR,1.0/3.0);
            return likelihoods;
        }
        return genotype.getLikelihoods().getAsMap(true);
    }
    //Returns the GenotypeType; returns UNVAILABLE if given null
    private GenotypeType getTypeSafeNull(Genotype genotype){
        if(genotype == null)
            return GenotypeType.UNAVAILABLE;
        return genotype.getType();
    }

    private static final double maxPL(double[] GLs) {
        double adjust = -1.0D / 0.0;
        double[] var3 = GLs;
        int var4 = GLs.length;

        for(int var5 = 0; var5 < var4; ++var5) {
            double l = var3[var5];
            adjust = Math.max(adjust, l);
        }

        return adjust;
    }
    public enum GenotypeTypeShort {
        HOM_REF,
        HET,
        HOM_VAR;

        private GenotypeTypeShort() {
        }
    }
    private int phaseTrioGenotypes(Allele ref, Allele alt, Genotype mother, Genotype father, Genotype child,ArrayList<Genotype> finalGenotypes, ReadLikelihoods<Allele> forMother, ReadLikelihoods<Haplotype> baseLevel, PloidyModel pl) {

        //Check whether it is  a pair or trio
        //Always assign the first parent as the parent having genotype information in pairs
        //Always assign the mother as the first parent in trios
        int parentsCalled = 0;
        Map<GenotypeType,Double> firstParentLikelihoods;
        Map<GenotypeType,Double> secondParentLikelihoods;
        EnumMap<GenotypeType, EnumMap<GenotypeType, EnumMap<GenotypeTypeShort, Double>>> abortusCummulative;

        abortusCummulative = new EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeTypeShort,Double>>>(GenotypeType.class);
        for(GenotypeType m : GenotypeType.values()){
            abortusCummulative.put(m,new EnumMap<GenotypeType,EnumMap<GenotypeTypeShort,Double>>(GenotypeType.class));
            for(GenotypeType f: GenotypeType.values()){
                abortusCummulative.get(m).put(f,new EnumMap<GenotypeTypeShort,Double>(GenotypeTypeShort.class));
                for(GenotypeTypeShort c : GenotypeTypeShort.values()){
                    abortusCummulative.get(m).get(f).put(c,0.0);
                }
            }
        }

        ArrayList<GenotypeType> bestFirstParentGenotype = new ArrayList<GenotypeType>();
        ArrayList<GenotypeType> bestSecondParentGenotype = new ArrayList<GenotypeType>();
        ArrayList<GenotypeType> bestChildGenotype = new ArrayList<GenotypeType>();
        GenotypeType pairSecondParentGenotype = null;
        System.out.println("mother " + mother.getAlleles());
        System.out.println("mother " + mother.getGenotypeString());
        System.out.println("father " + father.getAlleles());
        System.out.println("father " + father.getGenotypeString());

        if(mother == null || !mother.isCalled()){
            firstParentLikelihoods = getLikelihoodsAsMapSafeNull(father);
            secondParentLikelihoods = getLikelihoodsAsMapSafeNull(mother);
            bestFirstParentGenotype.add(getTypeSafeNull(father));
            bestSecondParentGenotype.add(getTypeSafeNull(mother));
            pairSecondParentGenotype = mother == null ? GenotypeType.UNAVAILABLE : mother.getType();
            if(father != null && father.isCalled())
                parentsCalled = 1;
        }
        else{

            firstParentLikelihoods = getLikelihoodsAsMapSafeNull(mother);
            secondParentLikelihoods = getLikelihoodsAsMapSafeNull(father);
            bestFirstParentGenotype.add(getTypeSafeNull(mother));
            bestSecondParentGenotype.add(getTypeSafeNull(father));
            if(father == null || !father.isCalled()){
                parentsCalled = 1;
                pairSecondParentGenotype = father == null ? GenotypeType.UNAVAILABLE : father.getType();
            }else{
                parentsCalled = 2;
            }
        }
        Map<GenotypeType,Double> childLikelihoods = getLikelihoodsAsMapSafeNull(child);
        bestChildGenotype.add(getTypeSafeNull(child));

        //Prior vars
        double bestConfigurationLikelihood = 0.0;
        double norm = 0.0;
        int configuration_index =0;
        ArrayList<Integer> bestMVCount = new ArrayList<Integer>();
        bestMVCount.add(0);

        GenotypingModel genotypingModel = new InfiniteRandomMatingPopulationModel();

        //Get the most likely combination
        //Only check for most likely combination if at least a parent and the child have genotypes
        if(child.isCalled() && parentsCalled > 0){
            int mvCount;
            int cumulativeMVCount = 0;
            double configurationLikelihood = 0;
            for(Map.Entry<GenotypeType,Double> childGenotype : childLikelihoods.entrySet()) {
                for (Map.Entry<GenotypeType, Double> firstParentGenotype : firstParentLikelihoods.entrySet()) {
                    for (Map.Entry<GenotypeType, Double> secondParentGenotype : secondParentLikelihoods.entrySet()) {
                        mvCount = mvCountMatrix.get(firstParentGenotype.getKey()).get(secondParentGenotype.getKey()).get(childGenotype.getKey());
                        /*TrioPhase tp = transmissionMatrix.get(secondParentGenotype.getKey()).get(firstParentGenotype.getKey()).get(childGenotype.getKey());*/

                        //System.out.println("transmissionMatrix = " + tp.);
                        //For parent/child pairs, sum over the possible genotype configurations of the missing parent
                        if (parentsCalled < 2) {
                            cumulativeMVCount += mvCount;
                            configurationLikelihood += mvCount > 0 ? Math.pow(deNovoPrior, mvCount) * firstParentGenotype.getValue() * secondParentGenotype.getValue() * childGenotype.getValue() : (1.0 - 11 * deNovoPrior) * firstParentGenotype.getValue() * secondParentGenotype.getValue() * childGenotype.getValue();
                        }
                        //Evaluate configurations of trios
                        else {
                            TrioPhase phasedTrioGenotypes;
                            phasedTrioGenotypes = transmissionMatrix.get(firstParentGenotype.getKey()).get(secondParentGenotype.getKey()).get(childGenotype.getKey());
                            System.out.println("phasedTrioGenotypes "
                                    + firstParentGenotype.getKey() + "," +
                                    secondParentGenotype.getKey() + "," +
                                    childGenotype.getKey() + "= ");

                            //double log10frequency_ploidy2 = 2;
                            //double log10_admixture_mother = 0.5;//0.2
                            //double log10_admixture_child = 0.5;//0.8

                            double log10frequency_ploidy2 = 0.301;
                            double log10_admixture_mother = -0.301;//0.2
                            double log10_admixture_child = -0.301;//0.8

                            //double log10_admixture_mother = -0.69897;//0.2
                            //double log10_admixture_child = -0.09691;//0.8

                            //double log10_admixture_mother = -1.301;
                            //double log10_admixture_child = -0.02227639;

                            //working with mother
                            double abortus_likehood_big = 0.0;
                            for (int readIdx = 0; readIdx < forMother.readCount(); readIdx++) {
                                System.out.println("readIdx = " + readIdx + " -- " + forMother.sampleReads(0).get(readIdx).getReadString());
                                double abortus_likehood = 0.0;
                                if (firstParentGenotype.getKey() == GenotypeType.HOM_REF ||
                                        firstParentGenotype.getKey() == GenotypeType.HOM_VAR) {
                                    Allele one_and_only = childGenotype.getKey() == GenotypeType.HOM_REF ? ref : alt;
                                    System.out.println("!!mother " + one_and_only + "" + one_and_only);


                                    for (int allele_idx = 0; allele_idx < forMother.alleles().size(); allele_idx++) {
                                        if (forMother.alleleAt(allele_idx).equals(one_and_only, true)) {
                                            System.out.println("+=" + forMother.valuesBySampleIndex[0][allele_idx][readIdx] +
                                                    " + " + log10frequency_ploidy2 + " + " + log10_admixture_mother);
                                            //abortus_likehood +=
                                            //            Math.pow(10, forMother.valuesBySampleIndex[0][allele_idx][readIdx])* log10frequency_ploidy2 * log10_admixture_mother;
                                            double lk_xAllels_admixture = (forMother.valuesBySampleIndex[0][allele_idx][readIdx] + log10frequency_ploidy2 + log10_admixture_mother);
                                            if (abortus_likehood == 0.0) {
                                                abortus_likehood = lk_xAllels_admixture;
                                            } else {
                                                abortus_likehood = MathUtils.approximateLog10SumLog10(abortus_likehood, lk_xAllels_admixture);
                                            }
                                        }

                                    }
                                } else if (firstParentGenotype.getKey() == GenotypeType.HET) {

                                    Allele mother_ref = ref;//mother.getAllele(0);
                                    Allele mother_var = alt;//mother.getAllele(1);
                                    System.out.println("!!mother " + mother_ref + "" + mother_var);
                                    int counter = 0;
                                    double l1 = 0.0, l2 = 0.0;
                                    for (int allele_idx = 0; allele_idx < forMother.alleles().size(); allele_idx++) {
                                        if (forMother.alleleAt(allele_idx).equals(mother_ref, true)) {
                                            counter++;
                                            l1 = forMother.valuesBySampleIndex[0][allele_idx][readIdx];
                                            System.out.println("l1=" + l1);
                                        }
                                        if (forMother.alleleAt(allele_idx).equals(mother_var, true)) {
                                            counter++;
                                            l2 = forMother.valuesBySampleIndex[0][allele_idx][readIdx];
                                            System.out.println("l2=" + l2);
                                        }
                                    }
                                    if (counter != 2) {
                                        System.out.println("PROBLEM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                                    } else {
                                        //abortus_likehood += Math.pow(10, MathUtils.approximateLog10SumLog10(l1,l2)) * log10_admixture_mother;
                                        double lk_xAllels_admixture = (MathUtils.approximateLog10SumLog10(l1, l2) + log10_admixture_mother);
                                        if (abortus_likehood == 0.0) {
                                            abortus_likehood = lk_xAllels_admixture;
                                        } else {
                                            abortus_likehood = MathUtils.approximateLog10SumLog10(abortus_likehood, lk_xAllels_admixture);
                                        }
                                    }
                                }

                                if (childGenotype.getKey() == GenotypeType.HOM_REF ||
                                        childGenotype.getKey() == GenotypeType.HOM_VAR) {
                                    //int combination_allele_idx = childGenotype.getKey() == GenotypeType.HOM_REF ? 0 : 1;
                                    Allele one_and_only = childGenotype.getKey() == GenotypeType.HOM_REF ? ref : alt;
                                    System.out.println("!!child " + one_and_only + "" + one_and_only);

                                    for (int allele_idx = 0; allele_idx < forMother.alleles().size(); allele_idx++) {
                                        if (forMother.alleleAt(allele_idx).equals(one_and_only, true)) {
                                            System.out.println("+=" + forMother.valuesBySampleIndex[0][allele_idx][readIdx] +
                                                    " + " + log10frequency_ploidy2 + " + " + log10_admixture_child);
                                            //abortus_likehood +=
                                            //        Math.pow(10, forMother.valuesBySampleIndex[0][allele_idx][readIdx])* log10frequency_ploidy2 * log10_admixture_child;
                                            double lk_xAllels_admixture = (forMother.valuesBySampleIndex[0][allele_idx][readIdx] + log10frequency_ploidy2 + log10_admixture_child);
                                            if (abortus_likehood == 0.0) {
                                                abortus_likehood = lk_xAllels_admixture;
                                            } else {
                                                abortus_likehood = MathUtils.approximateLog10SumLog10(abortus_likehood, lk_xAllels_admixture);
                                            }
                                        }

                                    }
                                } else if (childGenotype.getKey() == GenotypeType.HET) {
                                    Allele child_ref = ref;//child.getAllele(0);
                                    Allele child_var = alt;//child.getAllele(1);
                                    System.out.println("!!child " + child_ref + "" + child_var);
                                    int counter = 0;
                                    double l1 = 0.0, l2 = 0.0;
                                    for (int allele_idx = 0; allele_idx < forMother.alleles().size(); allele_idx++) {
                                        if (forMother.alleleAt(allele_idx).equals(child_ref, true)) {
                                            counter++;
                                            l1 = forMother.valuesBySampleIndex[0][allele_idx][readIdx];
                                            System.out.println("l1=" + l1);
                                        }
                                        if (forMother.alleleAt(allele_idx).equals(child_var, true)) {
                                            counter++;
                                            l2 = forMother.valuesBySampleIndex[0][allele_idx][readIdx];
                                            System.out.println("l2=" + l2);
                                        }
                                    }
                                    if (counter != 2) {
                                        System.out.println("PROBLEM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                                    } else {
                                        //abortus_likehood +=
                                        double lk_xAllels_admixture = (MathUtils.approximateLog10SumLog10(l1, l2) + log10_admixture_child);
                                        if (abortus_likehood == 0.0) {
                                            abortus_likehood = lk_xAllels_admixture;
                                        } else {
                                            abortus_likehood = MathUtils.approximateLog10SumLog10(abortus_likehood, lk_xAllels_admixture);
                                        }
                                        //abortus_likehood += Math.pow(10, MathUtils.approximateLog10SumLog10(l1,l2)) *  log10_admixture_child;
                                    }
                                }
                                abortus_likehood_big += (abortus_likehood - log10frequency_ploidy2);
                                /*if (abortus_likehood_big  == 0) {
                                    abortus_likehood_big = ;
                                } else {
                                    abortus_likehood_big *= abortus_likehood/log10frequency_ploidy2;
                                }*/
                                System.out.println("current abortus_likehood=" + abortus_likehood + ":abortus_likehood_big=" + abortus_likehood_big);
                            }
                            //forMother.valuesBySampleIndex[0][]
                            switch (childGenotype.getKey()) {
                                case HOM_REF:
                                    abortusCummulative.get(firstParentGenotype.getKey()).get(secondParentGenotype.getKey())
                                        .put(GenotypeTypeShort.HOM_REF, abortus_likehood_big);
                                    break;
                                case HET:
                                    abortusCummulative.get(firstParentGenotype.getKey()).get(secondParentGenotype.getKey())
                                            .put(GenotypeTypeShort.HET, abortus_likehood_big);
                                    break;
                                case HOM_VAR:
                                    abortusCummulative.get(firstParentGenotype.getKey()).get(secondParentGenotype.getKey())
                                            .put(GenotypeTypeShort.HOM_VAR, abortus_likehood_big);
                                    break;
                            }
                            ArrayList<Genotype> currentGenotypes = new ArrayList<>();
                            phasedTrioGenotypes.getPhasedGenotypes(ref, alt, mother, father, child, 0.0, currentGenotypes);
                            Genotype mother1 = currentGenotypes.get(0);
                            Genotype child1 = currentGenotypes.get(2);
                            System.out.println("Phasing " + ref.getBaseString() +
                                    " alt " + alt.getBaseString() +
                                    "mother" + mother +
                                    "father" + father +
                                    "child" + child + " -- "
                                    + child1.getAllele(1).getBaseString() + "/" + child1.getAllele(0).getBaseString() +
                                    ", " + mother1.getAllele(1).getBaseString() + "/" + mother1.getAllele(0).getBaseString());
                        }
                    }
                }
            }

            System.out.println("Normalizing abortus");
            for(GenotypeType m : GenotypeType.values()){
                for(GenotypeType f: GenotypeType.values()){
                    double[] gls = new double[3];
                    double[] pls = new double[3];
                    double[] final_p;

                    int pl_idx = 0;
                    for(GenotypeTypeShort c : GenotypeTypeShort.values()){
                        gls[pl_idx++] = abortusCummulative.get(m).get(f).get(c);

                    }
                    System.out.println("" + m + "," + f + "=" + gls);
                    double adjust = maxPL(gls);
                    for(int i = 0; i < gls.length; ++i) {
                        pls[i] = (int)Math.round(Math.min(-10.0D * (gls[i] - adjust), 2.147483647E9D));
                    }
                    System.out.println("" + m + "," + f + "=" + pls);
                    for(int i = 0; i < gls.length; ++i) {
                        pls[i] = pls[i]/(double)-10.0D;
                    }
                    final_p = GeneralUtils.normalizeFromLog10(pls);
                    System.out.println("" + m + "," + f + "=" + final_p);

                    pl_idx = 0;
                    for(GenotypeTypeShort c : GenotypeTypeShort.values()){
                        abortusCummulative.get(m).get(f).put(c, final_p[pl_idx++]);
                        System.out.println("" + m + "," + f + "," + c + "=" + final_p[pl_idx-1]);
                    }
                }
            }

            for(Map.Entry<GenotypeType,Double> childGenotype : childLikelihoods.entrySet()){
                for(Map.Entry<GenotypeType,Double> firstParentGenotype : firstParentLikelihoods.entrySet()) {
                    for (Map.Entry<GenotypeType, Double> secondParentGenotype : secondParentLikelihoods.entrySet()) {
                        mvCount = mvCountMatrix.get(firstParentGenotype.getKey()).get(secondParentGenotype.getKey()).get(childGenotype.getKey());
                        if (parentsCalled < 2) {
                            cumulativeMVCount += mvCount;
                            configurationLikelihood += mvCount > 0 ? Math.pow(deNovoPrior, mvCount) * firstParentGenotype.getValue() * secondParentGenotype.getValue() * childGenotype.getValue() : (1.0 - 11 * deNovoPrior) * firstParentGenotype.getValue() * secondParentGenotype.getValue() * childGenotype.getValue();
                        }
                        //Evaluate configurations of trios
                        else {
                            TrioPhase phasedTrioGenotypes;
                            phasedTrioGenotypes = transmissionMatrix.get(firstParentGenotype.getKey()).get(secondParentGenotype.getKey()).get(childGenotype.getKey());
                            //final GenotypingLikelihoods<Allele> likelihoods = genotypingModel.calculateLikelihoods(forMother, new GenotypingData<>(pl, forMother));
                          /*  GenotypingData gd = new GenotypingData<>(pl, forMother);
                            final AlleleListPermutation<Allele> permutation = AlleleListUtils.permutation(gd, forMother);
                            final AlleleLikelihoodMatrixMapper<Allele> alleleLikelihoodMatrixMapper = AlleleLikelihoodMatrixMapper.newInstance(permutation);

                            final int sampleCount = gd.sampleCount();

                            //singleSampleLikelihoods(forMother,gd,alleleLikelihoodMatrixMapper);

                            final PloidyModel ploidyModel = gd.ploidyModel();
                            final int samplePloidy = ploidyModel.samplePloidy(0);
                            final int alleleCount = forMother.alleleCount();
                            final GenotypeLikelihoodCalculator likelihoodsCalculator = GenotypeLikelihoodCalculators.getInstance(samplePloidy, alleleCount);
                            System.out.println("singleSampleLikelihoods (likelihoodsCalculator) = " + likelihoodsCalculator);
                            final ReadLikelihoods.Matrix<Allele> sampleLikelihoods = alleleLikelihoodMatrixMapper.map(gd.readLikelihoods().sampleMatrix(0));
                            System.out.println("singleSampleLikelihoods(sampleLikelihoods) = " + sampleLikelihoods);
                            //final List<GenotypeLikelihoods> genotypeLikelihoods = Collections.singletonList(likelihoodsCalculator.genotypeLikelihoods(sampleLikelihoods));
                            final int readCount = sampleLikelihoods.readCount();


                            likelihoodsCalculator.ensureReadCapacity(readCount);

                            /// [x][y][z] = z * LnLk(Read_x | Allele_y)
                            final double[] readLikelihoodComponentsByAlleleCount
                                    = likelihoodsCalculator.readLikelihoodComponentsByAlleleCount(sampleLikelihoods);
                            final double[][] genotypeLikelihoodByRead = likelihoodsCalculator.genotypeLikelihoodByRead(readLikelihoodComponentsByAlleleCount,readCount);
                            System.out.println("genotypeLikelihoodByRead.length=" + genotypeLikelihoodByRead.length);
                            for (int ij1 = 0; ij1 < genotypeLikelihoodByRead.length; ij1++ ) {
                                for (int ik1 = 0; ik1 < readCount; ik1++) {
                                    System.out.println(" genotypeLikelihoodByRead[" + ij1 + "] = " + genotypeLikelihoodByRead[ij1][ik1]);
                                }
                            }*/
                            //final double[] readLikelihoodsByGenotypeIndex = likelihoodsCalculator.genotypeLikelihoods(genotypeLikelihoodByRead, readCount);

                            //GenotypingLikelihoods<> likelihoods = new GenotypingLikelihoods<>(gd,ploidyModel,genotypeLikelihoods);
                            final double abortus_likehood_big;
                            GenotypeTypeShort gchild=GenotypeTypeShort.HOM_REF;
                                switch (childGenotype.getKey()) {
                                case HOM_REF:
                                    gchild = GenotypeTypeShort.HOM_REF;
                                    break;
                                case HET:
                                    gchild = GenotypeTypeShort.HET;
                                    break;
                                case HOM_VAR:
                                    gchild = GenotypeTypeShort.HOM_VAR;
                                    break;

                                    default: System.out.println("EXCEPTION " + childGenotype.getKey());
                                    break;
                            }
                            abortus_likehood_big = abortusCummulative.get(firstParentGenotype.getKey()).get(secondParentGenotype.getKey()).get(gchild);
                            System.out.println("Combo" + firstParentGenotype.getKey() + "," + secondParentGenotype.getKey() + "," +
                            gchild + " abortus_likehood_big=" + abortus_likehood_big);
                            //int iAllele = forMother.alleles().indexOf(child.getAllele(0));
                            //System.out.println("finded allele in forMother with index" + iAllele);
                            configurationLikelihood =  mvCount>0 ? Math.pow(deNovoPrior,mvCount)*firstParentGenotype.getValue()*secondParentGenotype.getValue()*abortus_likehood_big /**childGenotype.getValue()*/ : (1.0-11*deNovoPrior)*firstParentGenotype.getValue()*secondParentGenotype.getValue()*abortus_likehood_big/*childGenotype.getValue()*/;
                            System.out.println("likehoods = " + firstParentGenotype.getValue() +
                                    " second = " + secondParentGenotype.getValue() +
                                    "child = " + childGenotype.getValue());
                            //configurationLikelihood =  mvCount>0 ? Math.pow(deNovoPrior,mvCount)*firstParentGenotype.getValue()*secondParentGenotype.getValue()*childGenotype.getValue() : (1.0-11*deNovoPrior)*firstParentGenotype.getValue()*secondParentGenotype.getValue()*childGenotype.getValue();
                            norm += configurationLikelihood;
                            //Keep this combination if
                            //It has a better likelihood
                            //Or it has the same likelihood but requires less changes from original genotypes
                            if (configurationLikelihood > bestConfigurationLikelihood){
                                System.out.println("NEW BEST:" + firstParentGenotype + secondParentGenotype + childGenotype + " = " + configurationLikelihood + ">" + bestConfigurationLikelihood + " = " + bestFirstParentGenotype + bestSecondParentGenotype + bestChildGenotype );
                                bestConfigurationLikelihood = configurationLikelihood;
                                bestMVCount.clear();
                                bestMVCount.add(mvCount);
                                bestFirstParentGenotype.clear();
                                bestFirstParentGenotype.add(firstParentGenotype.getKey());
                                bestSecondParentGenotype.clear();
                                bestSecondParentGenotype.add(secondParentGenotype.getKey());
                                bestChildGenotype.clear();
                                bestChildGenotype.add(childGenotype.getKey());
                            }
                            else if(configurationLikelihood == bestConfigurationLikelihood) {
                                bestFirstParentGenotype.add(firstParentGenotype.getKey());
                                bestSecondParentGenotype.add(secondParentGenotype.getKey());
                                bestChildGenotype.add(childGenotype.getKey());
                                bestMVCount.add(mvCount);
                            }
                        }
                    }
                    //!!!CYCLE ENDS
                    //Evaluate configurations of parent/child pairs
                    if(parentsCalled<2){
                        norm += configurationLikelihood;
                        //Keep this combination if
                        //It has a better likelihood
                        //Or it has the same likelihood but requires less changes from original genotypes
                        if (configurationLikelihood > bestConfigurationLikelihood){
                            bestConfigurationLikelihood = configurationLikelihood;
                            bestMVCount.clear();
                            bestMVCount.add(cumulativeMVCount/3);
                            bestChildGenotype.clear();
                            bestFirstParentGenotype.clear();
                            bestSecondParentGenotype.clear();
                            bestChildGenotype.add(childGenotype.getKey());
                            bestFirstParentGenotype.add(firstParentGenotype.getKey());
                            bestSecondParentGenotype.add(pairSecondParentGenotype);
                        }
                        else if(configurationLikelihood == bestConfigurationLikelihood) {
                            bestFirstParentGenotype.add(firstParentGenotype.getKey());
                            bestSecondParentGenotype.add(pairSecondParentGenotype);
                            bestChildGenotype.add(childGenotype.getKey());
                            bestMVCount.add(cumulativeMVCount/3);
                        }
                        configurationLikelihood = 0;
                    }
                }
            }

            //normalize the best configuration probability
            bestConfigurationLikelihood = bestConfigurationLikelihood / norm;

            //In case of multiple equally likely combinations, take a random one
            if(bestFirstParentGenotype.size()>1){
                configuration_index = rand.nextInt(bestFirstParentGenotype.size()-1);
            }

        }
        else{
            bestConfigurationLikelihood = NO_TRANSMISSION_PROB;
        }

        TrioPhase phasedTrioGenotypes;
        if(parentsCalled < 2 && mother == null || !mother.isCalled())
            phasedTrioGenotypes = transmissionMatrix.get(bestSecondParentGenotype.get(configuration_index)).get(bestFirstParentGenotype.get(configuration_index)).get(bestChildGenotype.get(configuration_index));
        else
            phasedTrioGenotypes = transmissionMatrix.get(bestFirstParentGenotype.get(configuration_index)).get(bestSecondParentGenotype.get(configuration_index)).get(bestChildGenotype.get(configuration_index));

        //Return the phased genotypes
        phasedTrioGenotypes.getPhasedGenotypes(ref,alt,mother,father,child,bestConfigurationLikelihood,finalGenotypes);
        return bestMVCount.get(configuration_index);

    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    private final static List<VariantContext> NO_CALLS = Collections.emptyList();
    @Override
    public List<VariantContext> map( final ActiveRegion originalActiveRegion, final RefMetaDataTracker metaDataTracker ) {
        if ( justDetermineActiveRegions )
            // we're benchmarking ART and/or the active region determination code in the HC, just leave without doing any work
            return NO_CALLS;

        if (sampleNameToUse != null)
            removeReadsFromAllSamplesExcept(sampleNameToUse, originalActiveRegion);

        if( !originalActiveRegion.isActive() )
            // Not active so nothing to do!
            return referenceModelForNoVariation(originalActiveRegion, true);

        final List<VariantContext> givenAlleles = new ArrayList<>();
        if( HCAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES ) {
            for ( final VariantContext vc : metaDataTracker.getValues(HCAC.alleles) ) {
                if ( vc.isNotFiltered() ) {
                    givenAlleles.add(vc); // do something with these VCs during GGA mode
                }
            }
            // No alleles found in this region so nothing to do!
            if ( givenAlleles.isEmpty() ) { return referenceModelForNoVariation(originalActiveRegion, true); }
        } else {
            // No reads here so nothing to do!
            if( originalActiveRegion.size() == 0 ) { return referenceModelForNoVariation(originalActiveRegion, true); }
        }

        // run the local assembler, getting back a collection of information on how we should proceed
        final AssemblyResultSet untrimmedAssemblyResult = assembleReads(originalActiveRegion, givenAlleles);

        final TreeSet<VariantContext> allVariationEvents = untrimmedAssemblyResult.getVariationEvents();
        // TODO - line bellow might be unnecessary : it might be that assemblyResult will always have those alleles anyway
        // TODO - so check and remove if that is the case:
        allVariationEvents.addAll(givenAlleles);

        final ActiveRegionTrimmer.Result trimmingResult = trimmer.trim(originalActiveRegion,allVariationEvents);

        if (!trimmingResult.isVariationPresent() && !HCAC.disableOptimizations)
            return referenceModelForNoVariation(originalActiveRegion,false);

        final AssemblyResultSet assemblyResult =
                trimmingResult.needsTrimming() ? untrimmedAssemblyResult.trimTo(trimmingResult.getCallableRegion()) : untrimmedAssemblyResult;

        final ActiveRegion regionForGenotyping = assemblyResult.getRegionForGenotyping();

        if ( HCAC.bamWriter != null && HCAC.emitDroppedReads ) {
            haplotypeBAMWriter.addDroppedReadsFromDelta(DroppedReadsTracker.Reason.TRIMMMED, originalActiveRegion.getReads(), regionForGenotyping.getReads());
        }

        // filter out reads from genotyping which fail mapping quality based criteria
        //TODO - why don't do this before any assembly is done? Why not just once at the beginning of this method
        //TODO - on the originalActiveRegion?
        //TODO - if you move this up you might have to consider to change referenceModelForNoVariation
        //TODO - that does also filter reads.
        final Collection<GATKSAMRecord> filteredReads = filterNonPassingReads( regionForGenotyping );

        if ( HCAC.bamWriter != null && HCAC.emitDroppedReads ) {
            haplotypeBAMWriter.addDroppedReads(DroppedReadsTracker.Reason.FILTERED, filteredReads);
        }

        final Map<String, List<GATKSAMRecord>> perSampleFilteredReadList = splitReadsBySample( filteredReads );

        // abort early if something is out of the acceptable range
        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if( ! assemblyResult.isVariationPresent() && ! HCAC.disableOptimizations)
            return referenceModelForNoVariation(originalActiveRegion, false);

        // For sure this is not true if gVCF is on.
        if (dontGenotype) return NO_CALLS; // user requested we not proceed


        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if( regionForGenotyping.size() == 0 && ! HCAC.disableOptimizations) {
            // no reads remain after filtering so nothing else to do!
            return referenceModelForNoVariation(originalActiveRegion, false);
        }

        // evaluate each sample's reads against all haplotypes
        //logger.info("Computing read likelihoods with " + assemblyResult.regionForGenotyping.size() + " reads");
        final List<Haplotype> haplotypes = assemblyResult.getHaplotypeList();
        final Map<String,List<GATKSAMRecord>> reads = splitReadsBySample( regionForGenotyping.getReads() );

        // Calculate the likelihoods: CPU intensive part.
        final ReadLikelihoods<Haplotype> readLikelihoods =
                likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult,samplesList,reads);
        PloidyModel pl = new HomogeneousPloidyModel(samplesList, 2);
        for (int pi = regionForGenotyping.getLocation().getStart(); pi< regionForGenotyping.getLocation().getStop(); pi++) {
            System.out.println("finding ..." + pi);
            VariantContext vc = (VariantContext)metaDataTracker.getFirstValue(this.variantCollection.variants,
                    new GenomeLoc(regionForGenotyping.getLocation().getContig(),
                            regionForGenotyping.getLocation().getContigIndex(),
                            pi,
                            pi+1));
            if (vc == null)
            {
                continue;
                //return NO_CALLS;
            }
        //VariantContext vc = (VariantContext)metaDataTracker.getFirstValue(this.variantCollection.variants, regionForGenotyping.getLocation());
        //if (vc != null)
        //{
            System.out.println("vc = " + vc);
            System.out.println("Region for genotyping = " + regionForGenotyping.getLocation());
            Genotype mother = vc.getGenotype("mother");
            //Genotype abortus = vc.getGenotype("abortus05");
            Genotype abortus = vc.getGenotype("abortus");
            //Genotype abortus = vc.getGenotype("abortus02");
            Genotype father = vc.getGenotype("NA12892");
            System.out.println("" + mother + ", " + abortus + "," + father);
            System.out.println("hpls" + haplotypes);


            List<VariantContext> alls = new ArrayList();
            alls.add(vc);
            //AssemblyResultSet untrimmedForMother = assembleReads(originalActiveRegion, givenAlleles);
            AssemblyResultSet untrimmedForMother = null;
            try {
                untrimmedForMother = assembleReads(originalActiveRegion, alls);
            } catch (Exception e ) {
                System.out.println("Exception " + e.getLocalizedMessage());
                return NO_CALLS;
            }
           // System.out.println("alls = " + alls);

            TreeSet<VariantContext> motherVariationEvents = untrimmedForMother.getVariationEvents();
            System.out.println("untrimmedForMother =" + untrimmedForMother);
            System.out.println("motherVariationEvents(size) =" + motherVariationEvents.size());

            ActiveRegionTrimmer.Result trimmingResultMother = trimmer.trim(originalActiveRegion,motherVariationEvents);

            AssemblyResultSet assemblyMother =
                    trimmingResultMother.needsTrimming() ? untrimmedForMother.trimTo(trimmingResultMother.getCallableRegion()) : untrimmedForMother;


            /*priors test*/

            /*final VariantContextBuilder builder = new VariantContextBuilder(vc);
            VariantContext vc_familyPriors,vc_bothPriors;
            VariantContextUtils.calculateChromosomeCounts(builder, false);
            vc_familyPriors = builder.make();

            if (!skipPopulationPriors)
                vc_bothPriors = PosteriorLikelihoodsUtils.calculatePosteriorGLs(vc_familyPriors, otherVCs, missing * numRefIfMissing, globalPrior, !ignoreInputSamples, defaultToAC, useACoff);
            else {
                final VariantContextBuilder builder2 = new VariantContextBuilder(vc_familyPriors);
                VariantContextUtils.calculateChromosomeCounts(builder, false);
                vc_bothPriors = builder2.make();
            }
            vcfWriter.add(vc_bothPriors);*/

            /*priors test*/

            ReadLikelihoods<Haplotype> readLikelihoodsAbortusReadsForMother = this.likelihoodCalculationEngine.computeReadLikelihoods(assemblyMother, this.samplesList, reads);
            final List<Haplotype> haplotypes2 = assemblyMother.getHaplotypeList();
            //final Map<GATKSAMRecord,GATKSAMRecord> readRealignments2 = realignReadsToTheirBestHaplotype(readLikelihoodsAbortusReadsForMother, assemblyMother.getReferenceHaplotype(), assemblyMother.getPaddedReferenceLoc());
            /*for(int i =0; i<readLikelihoodsAbortusReadsForMother.alleleCount(); i++ ) {
                for(int j =0; j<readLikelihoodsAbortusReadsForMother.readCount(); j++ ) {
                    Allele a = readLikelihoodsAbortusReadsForMother.alleleAt(i);
                    System.out.println("Allele " + a.toString());

                    //for (VariantContext v : motherVariationEvents) {
                     //   for (Allele a1 : v.getAlleles()) {
                     //       if (a1.equals(a, true)) {
                    //            System.out.println("MATCH" + a + " == " + a1);
                    //        } else {
                    //            System.out.println("NO MATCH");
                    //        }
                    //    }
                    //}

                    System.out.println("--" + i + ":" + j + " = " + readLikelihoodsAbortusReadsForMother.valuesBySampleIndex[0][i][j]);
                }
            }*/
            // Realign reads to their best haplotype.
            //final Map<GATKSAMRecord,GATKSAMRecord> readRealignmentsForMother = realignReadsToTheirBestHaplotype(readLikelihoodsAbortusReadsForMother, assemblyMother.getReferenceHaplotype(), assemblyMother.getPaddedReferenceLoc());

            //readLikelihoodsAbortusReadsForMother.changeReads(readRealignmentsForMother);

            final Map<Integer, ReadLikelihoods<Allele>> calledHaplotypes2 =
            genotypingEngine.assignGenotypeLikelihoodsGiveMargs(
                    haplotypes2,
                    readLikelihoodsAbortusReadsForMother,
                    perSampleFilteredReadList,
                    assemblyMother.getFullReferenceWithPadding(),
                    assemblyMother.getPaddedReferenceLoc(),
                    regionForGenotyping.getLocation(),
                    getToolkit().getGenomeLocParser(),
                    metaDataTracker,
                    Collections.<VariantContext>emptyList(),
                    emitReferenceConfidence());
            /*HashMap<Allele, List<Haplotype>> oldToNew = new HashMap<>();

            ArrayList<VariantContext> varEventsList = new ArrayList<>();
            varEventsList.addAll(motherVariationEvents);

            final Map<HaplotypeCallerGenotypingEngine.Event, List<Haplotype>> eventMapper = HaplotypeCallerGenotypingEngine.createEventMapper(regionForGenotyping.getLocation().getLocation().getStart(), varEventsList, untrimmedForMother.getHaplotypeList());
            //make priority list
            System.out.println("createEventMapper");
            final List<String> priorityList = new LinkedList<>();
            for ( final VariantContext vcLs : varEventsList ) priorityList.add(vcLs.getSource());

            //make priority list
            System.out.println("priority list");
            VariantContext mergedVC = GATKVariantContextUtils.simpleMerge(varEventsList,
                    priorityList,
                    GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                    GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE,
                    false, false, null, false, false);
            System.out.println("mergedVC");
            final GenotypeLikelihoodsCalculationModel.Model calculationModel = mergedVC.isSNP() ? GenotypeLikelihoodsCalculationModel.Model.SNP
                    : GenotypeLikelihoodsCalculationModel.Model.INDEL;

            final Map<VariantContext, Allele> mergeMap = new LinkedHashMap<>();
            mergeMap.put(null, mergedVC.getReference()); // the reference event (null) --> the reference allele
            for(int iii = 0; iii < varEventsList.size(); iii++) {
                mergeMap.put(varEventsList.get(iii), mergedVC.getAlternateAllele(iii)); // BUGBUG: This is assuming that the order of alleles is the same as the priority list given to simpleMerge function
            }
            System.out.println("mergeMap");
            final Map<Allele, List<Haplotype>> alleleMapper = HaplotypeCallerGenotypingEngine.createAlleleMapper(mergeMap, eventMapper);

            System.out.println("alleleMapper");
            //mergedVC = removeAltAllelesIfTooManyGenotypes(ploidy, alleleMapper, mergedVC);

            //final ReadLikelihoods<Allele> readAlleleLikelihoods = readLikelihoods.marginalize(alleleMapper,
            //        genomeLocParser.createPaddedGenomeLoc(genomeLocParser.createGenomeLoc(mergedVC),
            //                ALLELE_EXTENSION));

            ReadLikelihoods<Allele> test = readLikelihoodsAbortusReadsForMother.marginalize(alleleMapper);
            System.out.println("marginalization " + test.readCount() + " -- " + oldToNew.values());

            for (Allele al : test.alleles()) {
                System.out.println("After marginalization = " + al.toString());
            }*/
            System.out.println("called2" + calledHaplotypes2);
            for (Map.Entry<Integer, ReadLikelihoods<Allele>>  i : calledHaplotypes2.entrySet()) {
                System.out.println("key = " + i.getKey() + " = " + i.getValue().alleles());
            }
            ArrayList s = (ArrayList)vc.getAttribute("PG");
            for (Object o: s) {
                System.out.println("PG(s):" + o) ;
            }
            System.out.println("Prior" + s.size() );
            for (VariantContext v : motherVariationEvents) {
                for (Allele a : v.getAlleles()) {
                    System.out.println("motherVariationEvents(i) = " + v + "---" + a.toString());
                }
            }
            System.out.println("readLikelihoodsAbortusReadsForMother(test) =" +
                    readLikelihoodsAbortusReadsForMother
                    + ":"
                    + readLikelihoodsAbortusReadsForMother.sampleCount()
                    + ":"
                    + readLikelihoodsAbortusReadsForMother.alleleCount()
                    + ":"
                    + readLikelihoodsAbortusReadsForMother.readCount());


            ArrayList<Genotype> trioGenotypes = new ArrayList<Genotype>(3);
            final int mvCount = phaseTrioGenotypes(vc.getReference(), vc.getAltAlleleWithHighestAlleleCount(), mother, father, abortus,trioGenotypes, calledHaplotypes2.get(pi), readLikelihoods, pl);
            System.out.println("Final:" + pi + " -- " + vc.getReference() + "," + vc.getAltAlleleWithHighestAlleleCount() + ", " +
                    mother + "," + father + "," + abortus + " = " + trioGenotypes);
            Genotype phasedMother = trioGenotypes.get(0);
            Genotype phasedFather = trioGenotypes.get(1);
            Genotype phasedChild = trioGenotypes.get(2);
            GenotypesContext genotypesContext = GenotypesContext.copy(vc.getGenotypes());
            genotypesContext.replace(phasedChild);
            if(mother != null){
                genotypesContext.replace(phasedMother);
                if(father != null){
                    genotypesContext.replace(phasedFather);
                }
            }
            else{
                genotypesContext.replace(phasedFather);
            }
            VariantContextBuilder builder = new VariantContextBuilder(vc);

            builder.genotypes(genotypesContext);
            vcfWriter.add(builder.make());
        }


        // Realign reads to their best haplotype.
        final Map<GATKSAMRecord,GATKSAMRecord> readRealignments = realignReadsToTheirBestHaplotype(readLikelihoods, assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc());

        if ( HCAC.bamWriter != null && HCAC.emitDroppedReads ) {
            haplotypeBAMWriter.addDroppedReadsFromDelta(
                    DroppedReadsTracker.Reason.REALIGNMENT_FAILURE,
                    regionForGenotyping.getReads(),
                    readRealignments.values());
        }

        readLikelihoods.changeReads(readRealignments);

        // Note: we used to subset down at this point to only the "best" haplotypes in all samples for genotyping, but there
        //  was a bad interaction between that selection and the marginalization that happens over each event when computing
        //  GLs.  In particular, for samples that are heterozygous non-reference (B/C) the marginalization for B treats the
        //  haplotype containing C as reference (and vice versa).  Now this is fine if all possible haplotypes are included
        //  in the genotyping, but we lose information if we select down to a few haplotypes.  [EB]

        final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes = genotypingEngine.assignGenotypeLikelihoods(
                haplotypes,
                readLikelihoods,
                perSampleFilteredReadList,
                assemblyResult.getFullReferenceWithPadding(),
                assemblyResult.getPaddedReferenceLoc(),
                regionForGenotyping.getLocation(),
                getToolkit().getGenomeLocParser(),
                metaDataTracker,
                (RTAC.consensusMode ? Collections.<VariantContext>emptyList() : givenAlleles),
                emitReferenceConfidence());
        for ( Haplotype h : calledHaplotypes.getCalledHaplotypes()) {
            System.out.println("called haplo = " + h.toString());
        }
        for ( VariantContext v1 : calledHaplotypes.getCalls()) {
            System.out.println("called calls = " + v1.toString());
        }
        System.out.println();
        if ( HCAC.bamWriter != null ) {
            final Set<Haplotype> calledHaplotypeSet = new HashSet<>(calledHaplotypes.getCalledHaplotypes());
            if (HCAC.disableOptimizations)
                calledHaplotypeSet.add(assemblyResult.getReferenceHaplotype());
            haplotypeBAMWriter.writeReadsAlignedToHaplotypes(
                    haplotypes,
                    assemblyResult.getPaddedReferenceLoc(),
                    haplotypes,
                    calledHaplotypeSet,
                    readLikelihoods);

            if ( HCAC.emitDroppedReads ) {
                haplotypeBAMWriter.writeDroppedReads();
            }
        }

        if( HCAC.DEBUG ) { logger.info("----------------------------------------------------------------------------------"); }


        if ( emitReferenceConfidence() ) {
            if ( !containsCalls(calledHaplotypes) ) {
                // no called all of the potential haplotypes
                return referenceModelForNoVariation(originalActiveRegion, false);
            } else {
                final List<VariantContext> result = new LinkedList<>();
                // output left-flanking non-variant section:
                if (trimmingResult.hasLeftFlankingRegion())
                    result.addAll(referenceModelForNoVariation(trimmingResult.nonVariantLeftFlankRegion(),false));
                // output variant containing region.
                result.addAll(referenceConfidenceModel.calculateRefConfidence(assemblyResult.getReferenceHaplotype(),
                        calledHaplotypes.getCalledHaplotypes(), assemblyResult.getPaddedReferenceLoc(), regionForGenotyping,
                        readLikelihoods, genotypingEngine.getPloidyModel(), genotypingEngine.getGenotypingModel(), calledHaplotypes.getCalls()));
                // output right-flanking non-variant section:
                if (trimmingResult.hasRightFlankingRegion())
                    result.addAll(referenceModelForNoVariation(trimmingResult.nonVariantRightFlankRegion(),false));
                return result;
            }
        } else
            return calledHaplotypes.getCalls();
    }

    /**
     * Returns a map with the original read as a key and the realigned read as the value.
     * <p>
     *     Missing keys or equivalent key and value pairs mean that the read was not realigned.
     * </p>
     * @return never {@code null}
     */
    private Map<GATKSAMRecord,GATKSAMRecord> realignReadsToTheirBestHaplotype(final ReadLikelihoods<Haplotype> originalReadLikelihoods, final Haplotype refHaplotype, final GenomeLoc paddedReferenceLoc) {

        final Collection<ReadLikelihoods<Haplotype>.BestAllele> bestAlleles = originalReadLikelihoods.bestAlleles();
        final Map<GATKSAMRecord,GATKSAMRecord> result = new HashMap<>(bestAlleles.size());

        for (final ReadLikelihoods<Haplotype>.BestAllele bestAllele : bestAlleles) {
            final GATKSAMRecord originalRead = bestAllele.read;
            final Haplotype bestHaplotype = bestAllele.allele;
            final boolean isInformative = bestAllele.isInformative();
            final GATKSAMRecord realignedRead = AlignmentUtils.createReadAlignedToRef(originalRead, bestHaplotype, refHaplotype, paddedReferenceLoc.getStart(), isInformative);
            result.put(originalRead,realignedRead);
        }
        return result;
    }

    private boolean containsCalls(final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes) {
        final List<VariantContext> calls = calledHaplotypes.getCalls();
        if (calls.isEmpty()) return false;
        for (final VariantContext call : calls)
            for (final Genotype genotype : call.getGenotypes())
                if (genotype.isCalled())
                    return true;
        return false;
    }

    /**
     * High-level function that runs the assembler on the active region reads,
     * returning a data structure with the resulting information needed
     * for further HC steps
     *
     * @param activeRegion the region we should assemble
     * @param giveAlleles additional alleles we might need to genotype (can be empty)
     * @return the AssemblyResult describing how to proceed with genotyping
     */
    protected AssemblyResultSet assembleReads(final ActiveRegion activeRegion, final List<VariantContext> giveAlleles) {
        // Create the reference haplotype which is the bases from the reference that make up the active region
        finalizeActiveRegion(activeRegion); // handle overlapping fragments, clip adapter and low qual tails
        if( HCAC.DEBUG ) { logger.info("Assembling " + activeRegion.getLocation() + " with " + activeRegion.size() + " reads:    (with overlap region = " + activeRegion.getExtendedLoc() + ")"); }

        final byte[] fullReferenceWithPadding = activeRegion.getActiveRegionReference(referenceReader, REFERENCE_PADDING);
        final GenomeLoc paddedReferenceLoc = getPaddedLoc(activeRegion);
        final Haplotype referenceHaplotype = createReferenceHaplotype(activeRegion, paddedReferenceLoc);

        // Create ReadErrorCorrector object if requested - will be used within assembly engine.
        ReadErrorCorrector readErrorCorrector = null;
        if (errorCorrectReads)
            readErrorCorrector = new ReadErrorCorrector(RTAC.kmerLengthForReadErrorCorrection, MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION, RTAC.minObservationsForKmerToBeSolid, HCAC.DEBUG, fullReferenceWithPadding);

        try {
            final AssemblyResultSet assemblyResultSet = assemblyEngine.runLocalAssembly( activeRegion, referenceHaplotype, fullReferenceWithPadding, paddedReferenceLoc, giveAlleles,readErrorCorrector );
            assemblyResultSet.debugDump(logger);
            return assemblyResultSet;

        } catch ( final Exception e ) {
            // Capture any exception that might be thrown, and write out the assembly failure BAM if requested
            if ( captureAssemblyFailureBAM ) {
                final SAMFileWriter writer = SAMFileWriterStub.createSAMFileWriter("assemblyFailure.bam", getToolkit());
                new DirectOutputTracker().addOutput((SAMFileWriterStub) writer);
                for ( final GATKSAMRecord read : activeRegion.getReads() ) {
                    writer.addAlignment(read);
                }
                writer.close();
            }
            throw e;
        }
    }

    /**
     * Helper function to create the reference haplotype out of the active region and a padded loc
     * @param activeRegion the active region from which to generate the reference haplotype
     * @param paddedReferenceLoc the GenomeLoc which includes padding and shows how big the reference haplotype should be
     * @return a non-null haplotype
     */
    private Haplotype createReferenceHaplotype(final ActiveRegion activeRegion, final GenomeLoc paddedReferenceLoc) {
        return ReferenceConfidenceModel.createReferenceHaplotype(activeRegion, activeRegion.getActiveRegionReference(referenceReader), paddedReferenceLoc);
    }

    /**
     * Create an ref model result (ref model or no calls depending on mode) for an active region without any variation
     * (not is active, or assembled to just ref)
     *
     * @param region the region to return a no-variation result
     * @param needsToBeFinalized should the region be finalized before computing the ref model (should be false if already done)
     * @return a list of variant contexts (can be empty) to emit for this ref region
     */
    private List<VariantContext> referenceModelForNoVariation(final ActiveRegion region, final boolean needsToBeFinalized) {
        if ( emitReferenceConfidence() ) {
            //TODO - why the activeRegion cannot manage its own one-time finalization and filtering?
            //TODO - perhaps we can remove the last parameter of this method and the three lines bellow?
            if ( needsToBeFinalized )
                finalizeActiveRegion(region);
            filterNonPassingReads(region);

            final GenomeLoc paddedLoc = region.getExtendedLoc();
            final Haplotype refHaplotype = createReferenceHaplotype(region, paddedLoc);
            final List<Haplotype> haplotypes = Collections.singletonList(refHaplotype);
            return referenceConfidenceModel.calculateRefConfidence(refHaplotype, haplotypes,
                    paddedLoc, region, createDummyStratifiedReadMap(refHaplotype, samplesList, region),
                    genotypingEngine.getPloidyModel(), genotypingEngine.getGenotypingModel(), Collections.<VariantContext>emptyList());
        } else
            return NO_CALLS;
    }

    /**
     * Create a context that maps each read to the reference haplotype with log10 L of 0
     * @param refHaplotype a non-null reference haplotype
     * @param samples a list of all samples
     * @param region the active region containing reads
     * @return a map from sample -> PerReadAlleleLikelihoodMap that maps each read to ref
     */
    public static ReadLikelihoods<Haplotype> createDummyStratifiedReadMap(final Haplotype refHaplotype,
                                                                          final SampleList samples,
                                                                          final ActiveRegion region) {
        return new ReadLikelihoods<>(samples, new IndexedAlleleList<>(refHaplotype),
                splitReadsBySample(samples, region.getReads()));
    }


    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(List<VariantContext> callsInRegion, Integer numCalledRegions) {
        /*for( final VariantContext call : callsInRegion ) {
            vcfWriter.add( call );
        }*/
        return (callsInRegion.isEmpty() ? 0 : 1) + numCalledRegions;
    }

    @Override
    public void onTraversalDone(Integer result) {
        genotypingEngine.printFinalMaxNumPLValuesWarning();
        if ( HCAC.emitReferenceConfidence == ReferenceConfidenceMode.GVCF ) ((GVCFWriter)vcfWriter).close(false); // GROSS -- engine forces us to close our own VCF writer since we wrapped it
        referenceConfidenceModel.close();
        //TODO remove the need to call close here for debugging, the likelihood output stream should be managed
        //TODO (open & close) at the walker, not the engine.
        likelihoodCalculationEngine.close();
        logger.info("Ran local assembly on " + result + " active regions");
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // private helper functions
    //
    //---------------------------------------------------------------------------------------------------------------

    private void finalizeActiveRegion( final ActiveRegion activeRegion ) {
        if (activeRegion.isFinalized()) return;

        // Loop through the reads hard clipping the adaptor and low quality tails
        final List<GATKSAMRecord> readsToUse = new ArrayList<>(activeRegion.getReads().size());
        for( final GATKSAMRecord myRead : activeRegion.getReads() ) {
            GATKSAMRecord clippedRead;
            if (errorCorrectReads)
                clippedRead = ReadClipper.hardClipLowQualEnds( myRead, MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION );
            else  // default case: clip low qual ends of reads
                clippedRead= ReadClipper.hardClipLowQualEnds( myRead, MIN_TAIL_QUALITY );

            if ( dontUseSoftClippedBases || ! ReadUtils.hasWellDefinedFragmentSize(clippedRead) ) {
                // remove soft clips if we cannot reliably clip off adapter sequence or if the user doesn't want to use soft clips at all
                clippedRead = ReadClipper.hardClipSoftClippedBases(clippedRead);
            } else {
                // revert soft clips so that we see the alignment start and end assuming the soft clips are all matches
                // TODO -- WARNING -- still possibility that unclipping the soft clips will introduce bases that aren't
                // TODO -- truly in the extended region, as the unclipped bases might actually include a deletion
                // TODO -- w.r.t. the reference.  What really needs to happen is that kmers that occur before the
                // TODO -- reference haplotype start must be removed
                clippedRead = ReadClipper.revertSoftClippedBases(clippedRead);
            }

            clippedRead = ( clippedRead.getReadUnmappedFlag() ? clippedRead : ReadClipper.hardClipAdaptorSequence( clippedRead ) );
            if( !clippedRead.isEmpty() && clippedRead.getCigar().getReadLength() > 0 ) {
                clippedRead = ReadClipper.hardClipToRegion( clippedRead, activeRegion.getExtendedLoc().getStart(), activeRegion.getExtendedLoc().getStop() );
                if( activeRegion.readOverlapsRegion(clippedRead) && clippedRead.getReadLength() > 0 ) {
                    //logger.info("Keeping read " + clippedRead + " start " + clippedRead.getAlignmentStart() + " end " + clippedRead.getAlignmentEnd());
                    readsToUse.add(clippedRead);
                }
            }
        }

        // TODO -- Performance optimization: we partition the reads by sample 4 times right now; let's unify that code.

        final List<GATKSAMRecord> downsampledReads = DownsamplingUtils.levelCoverageByPosition(ReadUtils.sortReadsByCoordinate(readsToUse), maxReadsInRegionPerSample, minReadsPerAlignmentStart);

        if ( HCAC.bamWriter != null && HCAC.emitDroppedReads ) {
            haplotypeBAMWriter.addDroppedReadsFromDelta(DroppedReadsTracker.Reason.DOWNSAMPLED, activeRegion.getReads(), downsampledReads);
        }

        // handle overlapping read pairs from the same fragment
        cleanOverlappingReadPairs(downsampledReads);

        activeRegion.clearReads();
        activeRegion.addAll(downsampledReads);
        activeRegion.setFinalized(true);
    }

    private Set<GATKSAMRecord> filterNonPassingReads( final ActiveRegion activeRegion ) {
        final Set<GATKSAMRecord> readsToRemove = new LinkedHashSet<>();
        for( final GATKSAMRecord rec : activeRegion.getReads() ) {
            if( rec.getReadLength() < READ_LENGTH_FILTER_THRESHOLD || rec.getMappingQuality() < READ_QUALITY_FILTER_THRESHOLD || (!isBadMateFilterDisabled && BadMateFilter.hasBadMate(rec)) || (keepRG != null && !rec.getReadGroup().getId().equals(keepRG)) ) {
                readsToRemove.add(rec);
            }
        }
        activeRegion.removeAll( readsToRemove );
        return readsToRemove;
    }

    private GenomeLoc getPaddedLoc( final ActiveRegion activeRegion ) {
        final int padLeft = Math.max(activeRegion.getExtendedLoc().getStart()-REFERENCE_PADDING, 1);
        final int padRight = Math.min(activeRegion.getExtendedLoc().getStop()+REFERENCE_PADDING, referenceReader.getSequenceDictionary().getSequence(activeRegion.getExtendedLoc().getContig()).getSequenceLength());
        return getToolkit().getGenomeLocParser().createGenomeLoc(activeRegion.getExtendedLoc().getContig(), padLeft, padRight);
    }

    private Map<String, List<GATKSAMRecord>> splitReadsBySample( final Collection<GATKSAMRecord> reads ) {
        return splitReadsBySample(samplesList, reads);
    }

    public static Map<String, List<GATKSAMRecord>> splitReadsBySample( final SampleList samplesList, final Collection<GATKSAMRecord> reads ) {
        final Map<String, List<GATKSAMRecord>> returnMap = new HashMap<>();
        final int sampleCount = samplesList.sampleCount();
        for (int i = 0; i < sampleCount; i++)
            returnMap.put(samplesList.sampleAt(i), new ArrayList<GATKSAMRecord>());

        for( final GATKSAMRecord read : reads )
            returnMap.get(read.getReadGroup().getSample()).add(read);

        return returnMap;
    }

    /**
     * Are we emitting a reference confidence in some form, or not?
     *
     * @return true if HC must emit reference confidence.
     */
    public boolean emitReferenceConfidence() {
        return HCAC.emitReferenceConfidence != ReferenceConfidenceMode.NONE;
    }

    /**
     * Clean up reads/bases that overlap within read pairs
     *
     * @param reads the list of reads to consider
     */
    private void cleanOverlappingReadPairs(final List<GATKSAMRecord> reads) {
        for ( final List<GATKSAMRecord> perSampleReadList : splitReadsBySample(reads).values() ) {
            final FragmentCollection<GATKSAMRecord> fragmentCollection = FragmentUtils.create(perSampleReadList);
            for ( final List<GATKSAMRecord> overlappingPair : fragmentCollection.getOverlappingPairs() )
                FragmentUtils.adjustQualsOfOverlappingPairedFragments(overlappingPair);
        }
    }

    private void removeReadsFromAllSamplesExcept(final String targetSample, final ActiveRegion activeRegion) {
        final Set<GATKSAMRecord> readsToRemove = new LinkedHashSet<>();
        for( final GATKSAMRecord rec : activeRegion.getReads() ) {
            if( !rec.getReadGroup().getSample().equals(targetSample) ) {
                readsToRemove.add(rec);
            }
        }
        activeRegion.removeAll( readsToRemove );

    }

    /**
     * Is writing to an output GVCF file?
     *
     * @return true if the VCF output file has a .g.vcf or .g.vcf.gz extension or if no output file
     */
    private boolean isGVCF() {
        final File file = ((VariantContextWriterStub) vcfWriter).getOutputFile();
        if ( file == null ){
            return true;
        } else {
            final String fileName = file.getName();
            return ( fileName.endsWith("." + GATKVCFUtils.GVCF_EXT) || fileName.endsWith("." + GATKVCFUtils.GVCF_GZ_EXT) );
        }
    }
}
