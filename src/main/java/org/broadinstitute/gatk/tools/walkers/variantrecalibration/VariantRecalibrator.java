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

package org.broadinstitute.gatk.tools.walkers.variantrecalibration;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.PartitionBy;
import org.broadinstitute.gatk.engine.walkers.PartitionType;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.R.RScriptExecutor;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.report.GATKReport;
import org.broadinstitute.gatk.utils.report.GATKReportColumn;
import org.broadinstitute.gatk.utils.report.GATKReportTable;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.gatk.utils.collections.ExpandingArrayList;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.io.Resource;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.*;

import Jama.Matrix;


import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.IOException;
/**
 * Build a recalibration model to score variant quality for filtering purposes
 *
 * <p>
 * The purpose of variant recalibration is to assign a well-calibrated probability to each variant call in a call set.
 * You can then create highly accurate call sets by filtering based on this single estimate for the accuracy of each call.
 * The approach taken by variant quality score recalibration is to develop a continuous, covarying estimate of the relationship
 * between SNP call annotations (such as QD, MQ, and ReadPosRankSum, for example) and the probability that a SNP is a true genetic
 * variant versus a sequencing or data processing artifact. This model is determined adaptively based on "true sites" provided
 * as input, typically HapMap 3 sites and those sites found to be polymorphic on the Omni 2.5M SNP chip array (in humans). This adaptive
 * error model can then be applied to both known and novel variation discovered in the call set of interest to evaluate the
 * probability that each call is real. The score that gets added to the INFO field of each variant is called the VQSLOD. It is
 * the log odds of being a true variant versus being false under the trained Gaussian mixture model.
 * </p>
 * 
 * <p>
 * This tool performs the first pass in a two-stage process called VQSR; the second pass is performed by the
 * <a href='https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_ApplyRecalibration.php'>ApplyRecalibration</a> tool.
 * In brief, the first pass consists of creating a Gaussian mixture model by looking at the distribution of annotation
 * values over a high quality subset of the input call set, and then scoring all input variants according to the model.
 * The second pass consists of filtering variants based on score cutoffs identified in the first pass.
 *</p>
 *
 * <p>VQSR is probably the hardest part of the Best Practices to get right, so be sure to read the
 * <a href='https://www.broadinstitute.org/gatk/guide/article?id=39'>method documentation</a>,
 * <a href='https://www.broadinstitute.org/gatk/guide/article?id=1259'>parameter recommendations</a> and
 * <a href='https://www.broadinstitute.org/gatk/guide/article?id=2805'>tutorial</a> to really understand what these
 * tools and how to use them for best results on your own data.</p>
 *
 * <h3>Inputs</h3>
 * <ul>
 * <li>The input raw variants to be recalibrated. These variant calls must be annotated with the annotations that will be used for modeling. 
 * If the calls come from multiple samples, they must have been obtained by joint calling the samples, either directly (running HaplotypeCaller 
 * on all samples together) or via the GVCF workflow (HaplotypeCaller with -ERC GVCF per-sample then GenotypeGVCFs on the resulting gVCFs) 
 * which is more scalable.</li>
 * <li>Known, truth, and training sets to be used by the algorithm. See the method documentation for more details.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 * <li>A recalibration table file that will be used by the ApplyRecalibration tool.</li>
 * <li>A tranches file which shows various metrics of the recalibration callset for slices of the data.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <p>Recalibrating SNPs in exome data:</p>
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -T VariantRecalibrator \
 *   -R reference.fasta \
 *   -input raw_variants.vcf \
 *   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.sites.vcf \
 *   -resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.b37.sites.vcf \
 *   -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.vcf
 *   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_135.b37.vcf \
 *   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff \
 *   -mode SNP \
 *   -recalFile output.recal \
 *   -tranchesFile output.tranches \
 *   -rscriptFile output.plots.R
 * </pre>
 *
 * <h3>Allele-specfic usage</h3>
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -T VariantRecalibrator \
 *   -R reference.fasta \
 *   -input raw_variants.withASannotations.vcf \
 *   -AS \
 *   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.sites.vcf \
 *   -resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.b37.sites.vcf \
 *   -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.vcf
 *   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_135.b37.vcf \
 *   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff \
 *   -mode SNP \
 *   -recalFile output.AS.recal \
 *   -tranchesFile output.AS.tranches \
 *   -rscriptFile output.plots.AS.R
 * </pre>
 * The input VCF must have been produced using allele-specific annotations in HaplotypeCaller.
 * Note that each allele will have a separate line in the output .recal file with its own VQSLOD and culprit that will be transferred to the final VCF in ApplyRecalibration.
 *
 * <h3>Caveats</h3>
 *
 * <ul>
 * <li>SNPs and indels must be recalibrated in separate runs (but it is not necessary to separate them into different files). Mixed records are treated as indels.</li>
 * <li>The values used in the example above are only meant to show how the command lines are composed.
 * They are not meant to be taken as specific recommendations of values to use in your own work, and they may be
 * different from the values cited elsewhere in our documentation. For the latest and greatest recommendations on
 * how to set parameter values for you own analyses, please read the Best Practices section of the documentation,
 * especially the <a href='https://www.broadinstitute.org/gatk/guide/article?id=1259'>FAQ document</a> on VQSR parameters.</li>
 * <li>Whole genomes and exomes take slightly different parameters, so make sure you adapt your commands accordingly! See the documents linked above for details.</li>
 * <li>If you work with small datasets (e.g. targeted capture experiments or small number of exomes), you will run into problems. Read the docs linked above for advice on how to deal with those issues.</li>
 * <li>In order to create the model reporting plots Rscript needs to be in your environment PATH (this is the scripting version of R, not the interactive version).
 * See <a target="r-project" href="http://www.r-project.org">http://www.r-project.org</a> for more info on how to download and install R.</li>
 * </ul>
 */

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARDISC, extraDocs = {CommandLineGATK.class} )
@PartitionBy(PartitionType.NONE)
public class VariantRecalibrator extends RodWalker<ExpandingArrayList<VariantDatum>, ExpandingArrayList<VariantDatum>> implements TreeReducible<ExpandingArrayList<VariantDatum>> {

    private static final String PLOT_TRANCHES_RSCRIPT = "plot_Tranches.R";

    @ArgumentCollection private VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();

    /////////////////////////////
    // Inputs
    /////////////////////////////
    /**
     * These variant calls must be annotated with the annotations that will be used for modeling. If the calls come from multiple samples, 
     * they must have been obtained by joint calling the samples, either directly (running HaplotypeCaller on all samples together) or 
     * via the GVCF workflow (HaplotypeCaller with -ERC GVCF per-sample then GenotypeGVCFs on the resulting gVCFs) which is more scalable. 
     * Note that the ability to pass multiple input files is only intended to facilitate scatter-gather parallelism (to enable e.g. running on VCFs 
     * generated per-chromosome), not to combine different callsets. The variant calls in the separate input files should not overlap. 
     */
    @Input(fullName="input", shortName = "input", doc="One or more VCFs of raw input variants to be recalibrated", required=true)
    public List<RodBindingCollection<VariantContext>> inputCollections;
    final private List<RodBinding<VariantContext>> input = new ArrayList<>();

    /**
     * These additional calls should be unfiltered and annotated with the error covariates that are intended to be used for modeling.
     */
    @Input(fullName="aggregate", shortName = "aggregate", doc="Additional raw input variants to be used in building the model", required=false)
    public List<RodBinding<VariantContext>> aggregate;

    /**
     * Any set of VCF files to use as lists of training, truth, or known sites.
     * Training - The program builds the Gaussian mixture model using input variants that overlap with these training sites.
     * Truth - The program uses these truth sites to determine where to set the cutoff in VQSLOD sensitivity.
     * Known - The program only uses known sites for reporting purposes (to indicate whether variants are already known or novel). They are not used in any calculations by the algorithm itself.
     * Bad - A database of known bad variants can be used to supplement the set of worst ranked variants (compared to the Gaussian mixture model) that the program selects from the data to model "bad" variants.
     */
    @Input(fullName="resource", shortName = "resource", doc="A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm (training and truth sets are required to run)", required=true)
    public List<RodBinding<VariantContext>> resource = Collections.emptyList();

    /////////////////////////////
    // Outputs
    /////////////////////////////
    @Output(fullName="recal_file", shortName="recalFile", doc="The output recal file used by ApplyRecalibration", required=true)
    protected VariantContextWriter recalWriter = null;

    @Output(fullName="tranches_file", shortName="tranchesFile", doc="The output tranches file used by ApplyRecalibration", required=true)
    protected File TRANCHES_FILE;

    /////////////////////////////
    // Additional Command Line Arguments
    /////////////////////////////
    /**
     * The expected transition / transversion ratio of true novel variants in your targeted region (whole genome, exome, specific
     * genes), which varies greatly by the CpG and GC content of the region. See expected Ti/Tv ratios section of the GATK best
     * practices documentation (http://www.broadinstitute.org/gatk/guide/best-practices) for more information.
     * Normal values are 2.15 for human whole genome values and 3.2 for human whole exomes. Note
     * that this parameter is used for display purposes only and isn't used anywhere in the algorithm!
     */
    @Argument(fullName="target_titv", shortName="titv", doc="The expected novel Ti/Tv ratio to use when calculating FDR tranches and for display on the optimization curve output figures. (approx 2.15 for whole genome experiments). ONLY USED FOR PLOTTING PURPOSES!", required=false)
    protected double TARGET_TITV = 2.15;

    /**
     * See the input VCF file's INFO field for a list of all available annotations.
     */
    @Argument(fullName="use_annotation", shortName="an", doc="The names of the annotations which should used for calculations", required=true)
    private List<String> USE_ANNOTATIONS = new ArrayList<String>();

    /**
     * Add truth sensitivity slices through the call set at the given values. The default values are 100.0, 99.9, 99.0, and 90.0
     * which will result in 4 estimated tranches in the final call set: the full set of calls (100% sensitivity at the accessible
     * sites in the truth set), a 99.9% truth sensitivity tranche, along with progressively smaller tranches at 99% and 90%.
     */
    @Argument(fullName="TStranche", shortName="tranche", doc="The levels of truth sensitivity at which to slice the data. (in percent, that is 1.0 for 1 percent)", required=false)
    private List<Double> TS_TRANCHES = new ArrayList<Double>(Arrays.asList(100.0, 99.9, 99.0, 90.0));
    /**
     * For this to work properly, the -ignoreFilter argument should also be applied to the ApplyRecalibration command.
     */
    @Argument(fullName="ignore_filter", shortName="ignoreFilter", doc="If specified, the variant recalibrator will also use variants marked as filtered by the specified filter name in the input VCF file", required=false)
    private List<String> IGNORE_INPUT_FILTERS = new ArrayList<String>();
    @Argument(fullName="ignore_all_filters", shortName="ignoreAllFilters", doc="If specified, the variant recalibrator will ignore all input filters. Useful to rerun the VQSR from a filtered output file.", required=false)
    private boolean IGNORE_ALL_FILTERS = false;
    @Output(fullName="rscript_file", shortName="rscriptFile", doc="The output rscript file generated by the VQSR to aid in visualization of the input data and learned model", required=false, defaultToStdout=false)
    private File RSCRIPT_FILE = null;

    /**
     *  This GATKReport gives information to describe the VQSR model fit. Normalized means for the positive model are
     *  concatenated as one table and negative model normalized means as another table. Covariances are also concatenated
     *  for positive and negative models, respectively. Tables of annotation means and standard deviations are provided
     *  to help describe the normalization. The model fit report can be read in with our R gsalib package. Individual
     *  model Gaussians can be subset by the value in the "Gaussian" column if desired.
     */
    @Argument(fullName="output_model", shortName = "outputModel", doc="If specified, the variant recalibrator will output the VQSR model to this file path.", required=false)
    private String outputModel = null;
    @Argument(fullName="input_model", shortName = "inputModel", doc="If specified, the variant recalibrator will read the VQSR model from this file path.", required=false)
    private String inputModel = "";
    //@Output(fullName="model_file", shortName = "modelFile", doc="A GATKReport containing the positive and negative model fits", required=false)
    //protected PrintStream modelReport = null;

    @Hidden
    @Argument(fullName="replicate", shortName="replicate", doc="Used to debug the random number generation inside the VQSR. Do not use.", required=false)
    protected int REPLICATE = 200;
    private ArrayList<Double> replicate = new ArrayList<>();

    /**
     * The statistical model being built by this tool may fail due to simple statistical sampling
     * issues. Rather than dying immediately when the initial model fails, this argument allows the
     * tool to restart with a different random seed and try to build the model again. The first
     * successfully built model will be kept.
     *
     * Note that the most common underlying cause of model building failure is that there is insufficient data to
     * build a really robust model. This argument provides a workaround for that issue but it is
     * preferable to provide this tool with more data (typically by including more samples or more territory)
     * in order to generate a more robust model.
     */
    @Advanced
    @Argument(fullName="max_attempts", shortName = "max_attempts", doc="Number of attempts to build a model before failing", required=false)
    protected int max_attempts = 1;

    /////////////////////////////
    // Debug Arguments
    /////////////////////////////
    @Advanced
    @Argument(fullName = "trustAllPolymorphic", shortName = "allPoly", doc = "Trust that all the input training sets' unfiltered records contain only polymorphic sites to drastically speed up the computation.", required = false)
    protected Boolean TRUST_ALL_POLYMORPHIC = false;

    @VisibleForTesting
    protected List<Integer> annotationOrder = null;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private VariantDataManager dataManager;
    private PrintStream tranchesStream;
    private final Set<String> ignoreInputFilterSet = new TreeSet<>();
    private final VariantRecalibratorEngine engine = new VariantRecalibratorEngine( VRAC );
    private GaussianMixtureModel goodModel = null;
    private GaussianMixtureModel badModel = null;

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public void initialize() {
        dataManager = new VariantDataManager( new ArrayList<>(USE_ANNOTATIONS), VRAC );

        if (RSCRIPT_FILE != null && !RScriptExecutor.RSCRIPT_EXISTS)
            Utils.warnUser(logger, String.format(
                    "Rscript not found in environment path. %s will be generated but PDF plots will not.",
                    RSCRIPT_FILE));

        if( IGNORE_INPUT_FILTERS != null ) {
            ignoreInputFilterSet.addAll( IGNORE_INPUT_FILTERS );
        }

        try {
            tranchesStream = new PrintStream(TRANCHES_FILE);
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(TRANCHES_FILE, e);
        }

        for( RodBinding<VariantContext> rod : resource ) {
            dataManager.addTrainingSet( new TrainingSet( rod ) );
        }

        if( !dataManager.checkHasTrainingSet() ) {
            throw new UserException.CommandLineException( "No training set found! Please provide sets of known polymorphic loci marked with the training=true ROD binding tag. For example, -resource:hapmap,VCF,known=false,training=true,truth=true,prior=12.0 hapmapFile.vcf" );
        }
        if( !dataManager.checkHasTruthSet() ) {
            throw new UserException.CommandLineException( "No truth set found! Please provide sets of known polymorphic loci marked with the truth=true ROD binding tag. For example, -resource:hapmap,VCF,known=false,training=true,truth=true,prior=12.0 hapmapFile.vcf" );
        }

        final File inputFile = new File(inputModel);
        if (inputFile.exists()) { // Load GMM from a file
            logger.info("Loading model from:" + inputModel);
            final GATKReport reportIn = new GATKReport(inputFile);

            // Read all the tables
            final GATKReportTable nmcTable = reportIn.getTable("NegativeModelCovariances");
            final GATKReportTable nmmTable = reportIn.getTable("NegativeModelMeans");
            final GATKReportTable nPMixTable = reportIn.getTable("BadGaussianPMix");
            final GATKReportTable pmcTable = reportIn.getTable("PositiveModelCovariances");
            final GATKReportTable pmmTable = reportIn.getTable("PositiveModelMeans");
            final GATKReportTable pPMixTable = reportIn.getTable("GoodGaussianPMix");
            final GATKReportTable anMeansTable = reportIn.getTable("AnnotationMeans");
            final GATKReportTable anStDevsTable = reportIn.getTable("AnnotationStdevs");

            orderAndValidateAnnotations(anMeansTable, dataManager.annotationKeys);

            final Map<String, Double> anMeans = getMapFromVectorTable(anMeansTable);
            final Map<String, Double> anStdDevs = getMapFromVectorTable(anStDevsTable);
            dataManager.setNormalization(anMeans, anStdDevs);

            goodModel = GMMFromTables(pmmTable, pmcTable, pPMixTable, annotationOrder.size());
            badModel = GMMFromTables(nmmTable, nmcTable, nPMixTable, annotationOrder.size());
        }

        final Set<VCFHeaderLine> hInfo = new HashSet<>();
        ApplyRecalibration.addVQSRStandardHeaderLines(hInfo);
        recalWriter.writeHeader( new VCFHeader(hInfo) );

        for( int iii = 0; iii < REPLICATE * 2; iii++ ) {
            replicate.add(Utils.getRandomGenerator().nextDouble());
        }

        // collect the actual rod bindings into a list for use later
        for ( final RodBindingCollection<VariantContext> inputCollection : inputCollections )
            input.addAll(inputCollection.getRodBindings());


    }

    /**
     * Order and validate annotations according to the annotations in the serialized model
     * Annotations on the command line must be the same as those in the model report or this will throw an exception.
     * Sets the {@code annotationOrder} list to map from command line order to the model report's order.
     * n^2 because we typically use 7 or less annotations.
     * @param annotationTable GATKReportTable of annotations read from the serialized model file
     */
    protected void orderAndValidateAnnotations(final GATKReportTable annotationTable, final List<String> annotationKeys){
        annotationOrder = new ArrayList<Integer>(annotationKeys.size());

        for (int i = 0; i < annotationTable.getNumRows(); i++){
            String serialAnno = (String)annotationTable.get(i, "Annotation");
            for (int j = 0; j < annotationKeys.size(); j++) {
                if (serialAnno.equals( annotationKeys.get(j) )){
                    annotationOrder.add(j);
                }
            }
        }

        if(annotationOrder.size() != annotationTable.getNumRows() || annotationOrder.size() != annotationKeys.size()) {
            final String errorMsg = "Annotations specified on the command line:"+annotationKeys.toString() +" do not match annotations in the model report:"+inputModel;
            throw new UserException.CommandLineException(errorMsg);
        }

    }


    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public ExpandingArrayList<VariantDatum> map( final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context ) {
        final ExpandingArrayList<VariantDatum> mapList = new ExpandingArrayList<>();

        if( tracker == null ) { // For some reason RodWalkers get map calls with null trackers
            return mapList;
        }

        mapList.addAll( addOverlappingVariants(input, true, tracker, context) );
        if( aggregate != null ) {
            mapList.addAll( addOverlappingVariants(aggregate, false, tracker, context) );
        }

        return mapList;
    }

    /**
     * Using the RefMetaDataTracker find overlapping variants and pull out the necessary information to create the VariantDatum
     * @param rods      the rods to search within
     * @param isInput   is this rod an -input rod?
     * @param tracker   the RefMetaDataTracker from the RODWalker map call
     * @param context   the AlignmentContext from the RODWalker map call
     * @return  a list of VariantDatums, can be empty
     */
    private List<VariantDatum> addOverlappingVariants( final List<RodBinding<VariantContext>> rods, final boolean isInput, final RefMetaDataTracker tracker, final AlignmentContext context ) {
        if( rods == null ) { throw new IllegalArgumentException("rods cannot be null."); }
        if( tracker == null ) { throw new IllegalArgumentException("tracker cannot be null."); }
        if( context == null ) { throw new IllegalArgumentException("context cannot be null."); }

        final ExpandingArrayList<VariantDatum> variants = new ExpandingArrayList<>();

        for( final VariantContext vc : tracker.getValues(rods, context.getLocation()) ) {
            if( vc != null && ( IGNORE_ALL_FILTERS || vc.isNotFiltered() || ignoreInputFilterSet.containsAll(vc.getFilters()) ) ) {
                if( VariantDataManager.checkVariationClass( vc, VRAC.MODE ) && !VRAC.useASannotations) {
                    addDatum(variants, isInput, tracker, context, vc, null, null);
                }
                else if(VRAC.useASannotations) {
                    for(final Allele allele : vc.getAlternateAlleles()) {
                        if(allele == Allele.SPAN_DEL)
                            continue;
                        if (VariantDataManager.checkVariationClass(vc, allele, VRAC.MODE))
                            addDatum(variants, isInput, tracker, context, vc, vc.getReference(), allele);
                    }
                }
            }
        }

        return variants;
    }

    /**
     * add a datum representing a variant site (or allele) to the data in {@code variants}, which represents the callset to be recalibrated
     * @param variants is modified by having a new VariantDatum added to it
     */
    private void addDatum(final ExpandingArrayList<VariantDatum> variants, final boolean isInput, final RefMetaDataTracker tracker, final AlignmentContext context, final VariantContext vc, final Allele refAllele, final Allele altAllele) {
        final VariantDatum datum = new VariantDatum();

        // Populate the datum with lots of fields from the VariantContext, unfortunately the VC is too big so we just pull in only the things we absolutely need.
        datum.referenceAllele = refAllele;
        datum.alternateAllele = altAllele;
        dataManager.decodeAnnotations(datum, vc, true); //BUGBUG: when run with HierarchicalMicroScheduler this is non-deterministic because order of calls depends on load of machine
        datum.loc = (isInput ? getToolkit().getGenomeLocParser().createGenomeLoc(vc) : null);
        datum.originalQual = vc.getPhredScaledQual();
        datum.isSNP = vc.isSNP() && vc.isBiallelic();
        datum.isTransition = datum.isSNP && GATKVariantContextUtils.isTransition(vc);
        datum.isAggregate = !isInput;

        // Loop through the training data sets and if they overlap this locus (and allele, if applicable) then update the prior and training status appropriately
        dataManager.parseTrainingSets(tracker, context.getLocation(), vc, datum, TRUST_ALL_POLYMORPHIC);
        final double priorFactor = QualityUtils.qualToProb(datum.prior);
        datum.prior = Math.log10(priorFactor) - Math.log10(1.0 - priorFactor);

        variants.add(datum);
    }


    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public ExpandingArrayList<VariantDatum> reduceInit() {
        return new ExpandingArrayList<>();
    }

    @Override
    public ExpandingArrayList<VariantDatum> reduce( final ExpandingArrayList<VariantDatum> mapValue, final ExpandingArrayList<VariantDatum> reduceSum ) {
        reduceSum.addAll( mapValue );
        return reduceSum;
    }

    @Override
    public ExpandingArrayList<VariantDatum> treeReduce( final ExpandingArrayList<VariantDatum> lhs, final ExpandingArrayList<VariantDatum> rhs ) {
        rhs.addAll( lhs );
        return rhs;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // on traversal done
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public void onTraversalDone( final ExpandingArrayList<VariantDatum> reduceSum ) {
        for (int i = 1; i <= max_attempts; i++) {
            try {
                dataManager.setData(reduceSum);
                dataManager.normalizeData(inputModel.isEmpty(), annotationOrder); // Each data point is now (x - mean) / standard deviation

                final List<VariantDatum> positiveTrainingData = dataManager.getTrainingData();
                final List<VariantDatum> negativeTrainingData;

                if (goodModel != null && badModel != null){ // GMMs were loaded from a file
                    logger.info("Using serialized GMMs from file...");
                    engine.evaluateData(dataManager.getData(), goodModel, false);
                    negativeTrainingData = dataManager.selectWorstVariants();
                } else { // Generate the GMMs from scratch
                    // Generate the positive model using the training data and evaluate each variant
                    goodModel = engine.generateModel(positiveTrainingData, VRAC.MAX_GAUSSIANS);
                    engine.evaluateData(dataManager.getData(), goodModel, false);
                    // Generate the negative model using the worst performing data and evaluate each variant contrastively
                    negativeTrainingData = dataManager.selectWorstVariants();
                    badModel = engine.generateModel(negativeTrainingData, Math.min(VRAC.MAX_GAUSSIANS_FOR_NEGATIVE_MODEL, VRAC.MAX_GAUSSIANS));

                    if (badModel.failedToConverge || goodModel.failedToConverge) {
                        throw new UserException("NaN LOD value assigned. Clustering with this few variants and these annotations is unsafe. Please consider " + (badModel.failedToConverge ? "raising the number of variants used to train the negative model (via --minNumBadVariants 5000, for example)." : "lowering the maximum number of Gaussians allowed for use in the model (via --maxGaussians 4, for example)."));
                    }

                }

                dataManager.dropAggregateData(); // Don't need the aggregate data anymore so let's free up the memory
                engine.evaluateData(dataManager.getData(), badModel, true);

                if (outputModel != null) {
                    try (PrintStream modelReporter = new PrintStream(outputModel)) {
                        GATKReport report = writeModelReport(goodModel, badModel, USE_ANNOTATIONS);
                        report.print(modelReporter);
                    } catch (FileNotFoundException e){
                        throw new UserException("Could not open output model file:" + outputModel);
                    }
                }

                engine.calculateWorstPerformingAnnotation(dataManager.getData(), goodModel, badModel);

                // Find the VQSLOD cutoff values which correspond to the various tranches of calls requested by the user
                final int nCallsAtTruth = TrancheManager.countCallsAtTruth(dataManager.getData(), Double.NEGATIVE_INFINITY);
                final TrancheManager.SelectionMetric metric = new TrancheManager.TruthSensitivityMetric(nCallsAtTruth);
                final List<Tranche> tranches = TrancheManager.findTranches(dataManager.getData(), TS_TRANCHES, metric, VRAC.MODE);
                tranchesStream.print(Tranche.tranchesString(tranches));

                logger.info("Writing out recalibration table...");
                dataManager.writeOutRecalibrationTable(recalWriter);
                if (RSCRIPT_FILE != null) {
                    logger.info("Writing out visualization Rscript file...");
                    createVisualizationScript(dataManager.getRandomDataForPlotting(1000, positiveTrainingData, negativeTrainingData, dataManager.getEvaluationData()), goodModel, badModel, 0.0, dataManager.getAnnotationKeys().toArray(new String[USE_ANNOTATIONS.size()]));
                }

                if (VRAC.MODE == VariantRecalibratorArgumentCollection.Mode.INDEL) {
                    // Print out an info message to make it clear why the tranches plot is not generated
                    logger.info("Tranches plot will not be generated since we are running in INDEL mode");
                } else {
                    // Execute the RScript command to plot the table of truth values
                    RScriptExecutor executor = new RScriptExecutor();
                    executor.addScript(new Resource(PLOT_TRANCHES_RSCRIPT, VariantRecalibrator.class));
                    executor.addArgs(TRANCHES_FILE.getAbsoluteFile(), TARGET_TITV);
                    // Print out the command line to make it clear to the user what is being executed and how one might modify it
                    logger.info("Executing: " + executor.getApproximateCommandLine());
                    executor.exec();
                }
                return;
            } catch (Exception e) {
                if (i == max_attempts) {
                    throw e;
                } else {
                    logger.info(String.format("Exception occurred on attempt %d of %d. Trying again. Message was: '%s'", i, max_attempts, e.getMessage()));
                }
            }
        }
    }

    /**
     * Rebuild a Gaussian Mixture Model from gaussian means and co-variates stored in a GATKReportTables
     * @param muTable           Table of Gaussian means
     * @param sigmaTable        Table of Gaussian co-variates
     * @param pmixTable         Table of PMixLog10 values
     * @param numAnnotations    Number of annotations, i.e. Dimension of the annotation space in which the Gaussians live
     * @return  a GaussianMixtureModel whose state reflects the state recorded in the tables.
     */
    protected GaussianMixtureModel GMMFromTables(final GATKReportTable muTable, final GATKReportTable sigmaTable, final GATKReportTable pmixTable, final int numAnnotations){
        List<MultivariateGaussian> gaussianList = new ArrayList<>();

        int curAnnotation = 0;
        for (GATKReportColumn reportColumn : muTable.getColumnInfo() ) {
            if (!reportColumn.getColumnName().equals("Gaussian")) {
                for (int row = 0; row < muTable.getNumRows(); row++) {
                    if (gaussianList.size() <= row){
                        MultivariateGaussian mg = new MultivariateGaussian(numAnnotations);
                        gaussianList.add(mg);
                    }
                    gaussianList.get(row).mu[curAnnotation] = (Double) muTable.get(row, reportColumn.getColumnName());
                }
                curAnnotation++;
            }
        }

        for (GATKReportColumn reportColumn : pmixTable.getColumnInfo() ) {
            if (reportColumn.getColumnName().equals("pMixLog10")) {
                for (int row = 0; row < pmixTable.getNumRows(); row++) {
                    gaussianList.get(row).pMixtureLog10 =  (Double) pmixTable.get(row, reportColumn.getColumnName());
                }
            }
        }

        int curJ = 0;
        for (GATKReportColumn reportColumn : sigmaTable.getColumnInfo() ) {
            if (reportColumn.getColumnName().equals("Gaussian")) continue;
            if (reportColumn.getColumnName().equals("Annotation")) continue;

            for (int row = 0; row < sigmaTable.getNumRows(); row++) {
                int curGaussian = row / numAnnotations;
                int curI = row % numAnnotations;
                double curVal = (Double) sigmaTable.get(row, reportColumn.getColumnName());
                gaussianList.get(curGaussian).sigma.set(curI, curJ, curVal);

            }
            curJ++;

        }

        return new GaussianMixtureModel(gaussianList, VRAC.SHRINKAGE, VRAC.DIRICHLET_PARAMETER, VRAC.PRIOR_COUNTS);

    }

    private Map<String, Double> getMapFromVectorTable(GATKReportTable vectorTable){
        Map<String, Double> dataMap = new HashMap<>();

        //do a row-major traversal
        for (int i = 0; i < vectorTable.getNumRows(); i++) {
            dataMap.put((String) vectorTable.get(i, 0), (Double) vectorTable.get(i, 1));
        }
        return dataMap;
    }

    protected GATKReport writeModelReport(final GaussianMixtureModel goodModel, final GaussianMixtureModel badModel, List<String> annotationList) {
        final String formatString = "%.16E";
        final GATKReport report = new GATKReport();

        if (dataManager != null) {  //for unit test
            final double[] meanVector = dataManager.getMeanVector();
            GATKReportTable annotationMeans = makeVectorTable("AnnotationMeans", "Mean for each annotation, used to normalize data", dataManager.annotationKeys, meanVector, "Mean", formatString);
            report.addTable(annotationMeans);

            final double[] varianceVector = dataManager.getVarianceVector();  //"varianceVector" is actually stdev
            GATKReportTable annotationVariances = makeVectorTable("AnnotationStdevs", "Standard deviation for each annotation, used to normalize data", dataManager.annotationKeys, varianceVector, "Standarddeviation", formatString);
            report.addTable(annotationVariances);
        }

        List<String> gaussianStrings = new ArrayList<>();
        final double[] pMixtureLog10s = new double[goodModel.getModelGaussians().size()];
        int idx = 0;

        for( final MultivariateGaussian gaussian : goodModel.getModelGaussians() ) {
            pMixtureLog10s[idx] = gaussian.pMixtureLog10;
            gaussianStrings.add(Integer.toString(idx++) );
        }

        GATKReportTable goodPMix = makeVectorTable("GoodGaussianPMix", "Pmixture log 10 used to evaluate model", gaussianStrings, pMixtureLog10s, "pMixLog10", formatString, "Gaussian");
        report.addTable(goodPMix);

        gaussianStrings.clear();
        final double[] pMixtureLog10sBad = new double[badModel.getModelGaussians().size()];
        idx = 0;

        for( final MultivariateGaussian gaussian : badModel.getModelGaussians() ) {
            pMixtureLog10sBad[idx] = gaussian.pMixtureLog10;
            gaussianStrings.add(Integer.toString(idx++));
        }
        GATKReportTable badPMix = makeVectorTable("BadGaussianPMix", "Pmixture log 10 used to evaluate model", gaussianStrings, pMixtureLog10sBad, "pMixLog10", formatString, "Gaussian");
        report.addTable(badPMix);


        //The model and Gaussians don't know what the annotations are, so get them from this class
        //VariantDataManager keeps the annotation in the same order as the argument list
        GATKReportTable positiveMeans = makeMeansTable("PositiveModelMeans", "Vector of annotation values to describe the (normalized) mean for each Gaussian in the positive model", annotationList, goodModel, formatString);
        report.addTable(positiveMeans);

        GATKReportTable positiveCovariance = makeCovariancesTable("PositiveModelCovariances", "Matrix to describe the (normalized) covariance for each Gaussian in the positive model; covariance matrices are joined by row", annotationList, goodModel, formatString);
        report.addTable(positiveCovariance);

        //do the same for the negative model means
        GATKReportTable negativeMeans = makeMeansTable("NegativeModelMeans", "Vector of annotation values to describe the (normalized) mean for each Gaussian in the negative model", annotationList, badModel, formatString);
        report.addTable(negativeMeans);

        GATKReportTable negativeCovariance = makeCovariancesTable("NegativeModelCovariances", "Matrix to describe the (normalized) covariance for each Gaussian in the negative model; covariance matrices are joined by row", annotationList, badModel, formatString);
        report.addTable(negativeCovariance);

        return report;
    }

    protected GATKReportTable makeVectorTable(final String tableName, final String tableDescription, final List<String> annotationList, final double[] perAnnotationValues, final String columnName, final String formatString) {
        return makeVectorTable(tableName, tableDescription, annotationList, perAnnotationValues, columnName, formatString, "Annotation");
    }

    protected GATKReportTable makeVectorTable(final String tableName, final String tableDescription, final List<String> annotationList, final double[] perAnnotationValues, final String columnName, final String formatString, final String firstColumn) {
        GATKReportTable vectorTable = new GATKReportTable(tableName, tableDescription, annotationList.size(), GATKReportTable.TableSortingWay.DO_NOT_SORT);
        vectorTable.addColumn(firstColumn);
        vectorTable.addColumn(columnName, formatString);
        for (int i = 0; i < perAnnotationValues.length; i++) {
            vectorTable.addRowIDMapping(annotationList.get(i), i, true);
            vectorTable.set(i, 1, perAnnotationValues[i]);
        }
        return vectorTable;
    }

    private GATKReportTable makeMeansTable(final String tableName, final String tableDescription, final List<String> annotationList, final GaussianMixtureModel model, final String formatString) {
        GATKReportTable meansTable = new GATKReportTable(tableName, tableDescription, annotationList.size(), GATKReportTable.TableSortingWay.DO_NOT_SORT);
        meansTable.addColumn("Gaussian");
        for (final String annotationName : annotationList) {
            meansTable.addColumn(annotationName, formatString);
        }
        final List<MultivariateGaussian> modelGaussians = model.getModelGaussians();
        for (int i = 0; i < modelGaussians.size(); i++) {
            final MultivariateGaussian gaussian = modelGaussians.get(i);
            final double[] meanVec = gaussian.mu;
            if (meanVec.length != annotationList.size())
                throw new IllegalStateException("Gaussian mean vector does not have the same size as the list of annotations");
            meansTable.addRowIDMapping(i, i, true);
            for (int j = 0; j < annotationList.size(); j++)
                meansTable.set(i, annotationList.get(j), meanVec[j]);
        }
        return meansTable;
    }

    private GATKReportTable makeCovariancesTable(final String tableName, final String tableDescription, final List<String> annotationList, final GaussianMixtureModel model, final String formatString) {
        GATKReportTable modelCovariances = new GATKReportTable(tableName, tableDescription, annotationList.size()+2, GATKReportTable.TableSortingWay.DO_NOT_SORT); //+2 is for Gaussian and Annotation columns
        modelCovariances.addColumn("Gaussian");
        modelCovariances.addColumn("Annotation");
        for (final String annotationName : annotationList) {
            modelCovariances.addColumn(annotationName, formatString);
        }
        final List<MultivariateGaussian> modelGaussians = model.getModelGaussians();
        for (int i = 0; i < modelGaussians.size(); i++) {
            final MultivariateGaussian gaussian = modelGaussians.get(i);
            final Matrix covMat = gaussian.sigma;
            if (covMat.getRowDimension() != annotationList.size() || covMat.getColumnDimension() != annotationList.size())
                throw new IllegalStateException("Gaussian covariance matrix does not have the same size as the list of annotations");
            for (int j = 0; j < annotationList.size(); j++) {
                modelCovariances.set(j + i * annotationList.size(), "Gaussian", i);
                modelCovariances.set(j + i * annotationList.size(), "Annotation", annotationList.get(j));
                for (int k = 0; k < annotationList.size(); k++) {
                    modelCovariances.set(j + i * annotationList.size(), annotationList.get(k), covMat.get(j, k));

                }
            }
        }
        return modelCovariances;
    }

    private void createVisualizationScript( final List<VariantDatum> randomData, final GaussianMixtureModel goodModel, final GaussianMixtureModel badModel, final double lodCutoff, final String[] annotationKeys ) {
        PrintStream stream;
        try {
            stream = new PrintStream(RSCRIPT_FILE);
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(RSCRIPT_FILE, e);
        }

        // We make extensive use of the ggplot2 R library: http://had.co.nz/ggplot2/
        stream.println("library(ggplot2)");
        // For compactPDF in R 2.13+
        stream.println("library(tools)");
        // For graphical functions R 2.14.2+
        stream.println("library(grid)");

        createArrangeFunction( stream );

        stream.println("outputPDF <- \"" + RSCRIPT_FILE + ".pdf\"");
        stream.println("pdf(outputPDF)"); // Unfortunately this is a huge pdf file, BUGBUG: need to work on reducing the file size

        for(int iii = 0; iii < annotationKeys.length; iii++) {
            for( int jjj = iii + 1; jjj < annotationKeys.length; jjj++) {
                logger.info( "Building " + annotationKeys[iii] + " x " + annotationKeys[jjj] + " plot...");

                final List<VariantDatum> fakeData = new ExpandingArrayList<>();
                double minAnn1 = 100.0, maxAnn1 = -100.0, minAnn2 = 100.0, maxAnn2 = -100.0;
                for( final VariantDatum datum : randomData ) {
                    minAnn1 = Math.min(minAnn1, datum.annotations[iii]);
                    maxAnn1 = Math.max(maxAnn1, datum.annotations[iii]);
                    minAnn2 = Math.min(minAnn2, datum.annotations[jjj]);
                    maxAnn2 = Math.max(maxAnn2, datum.annotations[jjj]);
                }
                // Create a fake set of data which spans the full extent of these two annotation dimensions in order to calculate the model PDF projected to 2D
                final double NUM_STEPS = 60.0;
                for(double ann1 = minAnn1; ann1 <= maxAnn1; ann1+= (maxAnn1 - minAnn1) / NUM_STEPS) {
                    for(double ann2 = minAnn2; ann2 <= maxAnn2; ann2+= (maxAnn2 - minAnn2) / NUM_STEPS) {
                        final VariantDatum datum = new VariantDatum();
                        datum.prior = 0.0;
                        datum.annotations = new double[randomData.get(0).annotations.length];
                        datum.isNull = new boolean[randomData.get(0).annotations.length];
                        for(int ann=0; ann< datum.annotations.length; ann++) {
                            datum.annotations[ann] = 0.0;
                            datum.isNull[ann] = true;
                        }
                        datum.annotations[iii] = ann1;
                        datum.annotations[jjj] = ann2;
                        datum.isNull[iii] = false;
                        datum.isNull[jjj] = false;
                        fakeData.add(datum);
                    }
                }

                engine.evaluateData( fakeData, goodModel, false );
                engine.evaluateData( fakeData, badModel, true );

                stream.print("surface <- c(");
                for( final VariantDatum datum : fakeData ) {
                    stream.print(String.format("%.4f, %.4f, %.4f, ",
                            dataManager.denormalizeDatum(datum.annotations[iii], iii),
                            dataManager.denormalizeDatum(datum.annotations[jjj], jjj),
                            Math.min(4.0, Math.max(-4.0, datum.lod))));
                }
                stream.println("NA,NA,NA)");
                stream.println("s <- matrix(surface,ncol=3,byrow=T)");

                stream.print("data <- c(");
                for( final VariantDatum datum : randomData ) {
                    stream.print(String.format("%.4f, %.4f, %.4f, %d, %d,",
                            dataManager.denormalizeDatum(datum.annotations[iii], iii),
                            dataManager.denormalizeDatum(datum.annotations[jjj], jjj),
                            (datum.lod < lodCutoff ? -1.0 : 1.0),
                            (datum.atAntiTrainingSite ? -1 : (datum.atTrainingSite ? 1 : 0)), (datum.isKnown ? 1 : -1)));
                }
                stream.println("NA,NA,NA,NA,1)");
                stream.println("d <- matrix(data,ncol=5,byrow=T)");

                final String surfaceFrame = "sf." + annotationKeys[iii] + "." + annotationKeys[jjj];
                final String dataFrame = "df." + annotationKeys[iii] + "." + annotationKeys[jjj];

                stream.println(surfaceFrame + " <- data.frame(x=s[,1], y=s[,2], lod=s[,3])");
                stream.println(dataFrame + " <- data.frame(x=d[,1], y=d[,2], retained=d[,3], training=d[,4], novelty=d[,5])");
                stream.println("dummyData <- " + dataFrame + "[1,]");
                stream.println("dummyData$x <- NaN");
                stream.println("dummyData$y <- NaN");
                stream.println("p <- ggplot(data=" + surfaceFrame + ", aes(x=x, y=y)) +theme(panel.background = element_rect(fill = \"white\"), panel.grid.minor = element_line(colour = \"white\"), panel.grid.major = element_line(colour = \"white\"))");
                stream.println("p1 = p +ggtitle(\"model PDF\") + labs(x=\""+ annotationKeys[iii] +"\", y=\""+ annotationKeys[jjj] +"\") + geom_tile(aes(fill = lod)) + scale_fill_gradient(high=\"green\", low=\"red\", space=\"rgb\")");
                stream.println("p <- qplot(x,y,data=" + dataFrame + ", color=retained, alpha=I(1/7),legend=FALSE) +theme(panel.background = element_rect(fill = \"white\"), panel.grid.minor = element_line(colour = \"white\"), panel.grid.major = element_line(colour = \"white\"))");
                stream.println("q <- geom_point(aes(x=x,y=y,color=retained),data=dummyData, alpha=1.0, na.rm=TRUE)");
                stream.println("p2 = p + q + labs(x=\""+ annotationKeys[iii] +"\", y=\""+ annotationKeys[jjj] +"\") + scale_colour_gradient(name=\"outcome\", high=\"black\", low=\"red\",breaks=c(-1,1),guide=\"legend\",labels=c(\"filtered\",\"retained\"))");
                stream.println("p <- qplot(x,y,data="+ dataFrame + "["+dataFrame+"$training != 0,], color=training, alpha=I(1/7)) +theme(panel.background = element_rect(fill = \"white\"), panel.grid.minor = element_line(colour = \"white\"), panel.grid.major = element_line(colour = \"white\"))");
                stream.println("q <- geom_point(aes(x=x,y=y,color=training),data=dummyData, alpha=1.0, na.rm=TRUE)");
                stream.println("p3 = p + q + labs(x=\""+ annotationKeys[iii] +"\", y=\""+ annotationKeys[jjj] +"\") + scale_colour_gradient(high=\"green\", low=\"purple\",breaks=c(-1,1),guide=\"legend\", labels=c(\"neg\", \"pos\"))");
                stream.println("p <- qplot(x,y,data=" + dataFrame + ", color=novelty, alpha=I(1/7)) +theme(panel.background = element_rect(fill = \"white\"), panel.grid.minor = element_line(colour = \"white\"), panel.grid.major = element_line(colour = \"white\"))");
                stream.println("q <- geom_point(aes(x=x,y=y,color=novelty),data=dummyData, alpha=1.0, na.rm=TRUE)");
                stream.println("p4 = p + q + labs(x=\""+ annotationKeys[iii] +"\", y=\""+ annotationKeys[jjj] +"\") + scale_colour_gradient(name=\"novelty\", high=\"blue\", low=\"red\",breaks=c(-1,1),guide=\"legend\", labels=c(\"novel\",\"known\"))");
                stream.println("arrange(p1, p2, p3, p4, ncol=2)");
            }
        }
        stream.println("dev.off()");

        stream.println("if (exists(\"compactPDF\")) {");
        stream.println("compactPDF(outputPDF)");
        stream.println("}");

        stream.close();

        // Execute Rscript command to generate the clustering plots
        RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(RSCRIPT_FILE);
        logger.info("Executing: " + executor.getApproximateCommandLine());
        executor.exec();
     }

    // The Arrange function is how we place the 4 model plots on one page
    // from http://gettinggeneticsdone.blogspot.com/2010/03/arrange-multiple-ggplot2-plots-in-same.html
    private void createArrangeFunction( final PrintStream stream ) {
        stream.println("vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)");
        stream.println("arrange <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {");
        stream.println("dots <- list(...)");
        stream.println("n <- length(dots)");
        stream.println("if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}");
        stream.println("if(is.null(nrow)) { nrow = ceiling(n/ncol)}");
        stream.println("if(is.null(ncol)) { ncol = ceiling(n/nrow)}");
        stream.println("grid.newpage()");
        stream.println("pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )");
        stream.println("ii.p <- 1");
        stream.println("for(ii.row in seq(1, nrow)){");
        stream.println("ii.table.row <- ii.row ");
        stream.println("if(as.table) {ii.table.row <- nrow - ii.table.row + 1}");
        stream.println("for(ii.col in seq(1, ncol)){");
        stream.println("ii.table <- ii.p");
        stream.println("if(ii.p > n) break");
        stream.println("print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))");
        stream.println("ii.p <- ii.p + 1");
        stream.println("}");
        stream.println("}");
        stream.println("}");
    }
}
