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

import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.io.FileUtils;
import org.apache.log4j.Level;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.variant.GATKVCFIndexType;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

public class HaplotypeCallerGVCFIntegrationTest extends WalkerTest {

    @DataProvider(name = "MyDataProviderHaploid")
    public Object[][] makeMyDataProviderHaploid() {
        List<Object[]> tests = new ArrayList<>();

        final String PCRFreeIntervals = "-L 20:10,000,000-10,010,000";
        final String WExIntervals = "-L 20:10,000,000-10,100,000 -isr INTERSECTION -L " + hg19Chr20Intervals;

        // this functionality can be adapted to provide input data for whatever you might want in your data
        //TODO the latest independent exact AC calculation for haploids causes a clear variant to be lost here.
        //TODO this might need to be addressed at some point.
        //TODO the following test is commented out for the record
        //tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.NONE, PCRFreeIntervals, "7f09c261950bf86e435edfa69ed2ec71"});
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.NONE, PCRFreeIntervals, "db25572025ab2481b40f100ec499db45"});
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.BP_RESOLUTION, PCRFreeIntervals, "568ae19bd27bf81add12b8fd53a37b30"});
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.GVCF, PCRFreeIntervals, "57ac809fc00eb87724eef6064193b2c9"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.NONE, WExIntervals, "2870f9cf471a737a5d6b23ba8e250e03"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.BP_RESOLUTION, WExIntervals, "516e51a94e968cbffdc37f6da470677c"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.GVCF, WExIntervals, "e6c8108166c5c112e29bf8310e141dfe"});

        return tests.toArray(new Object[][]{});
    }


    @DataProvider(name = "MyDataProvider")
    public Object[][] makeMyDataProvider() {
        List<Object[]> tests = new ArrayList<>();

        final String PCRFreeIntervals = "-L 20:10,000,000-10,010,000";
        final String WExIntervals = "-L 20:10,000,000-10,100,000 -isr INTERSECTION -L " + hg19Chr20Intervals;

        // this functionality can be adapted to provide input data for whatever you might want in your data
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.NONE, PCRFreeIntervals, "77c2e2075895ea7635949534d38b45c5"});
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.BP_RESOLUTION, PCRFreeIntervals, "5f8be6e9bdf2e8c06ed31b71804c8edc"});
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.GVCF, PCRFreeIntervals, "7a368e838fbb17ccfa529ba9210b02cd"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.NONE, WExIntervals, "a08926513178b6f82c37bea85e847437"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.BP_RESOLUTION, WExIntervals, "7debb053c4241c070625d7041da71196"});

        final String NA12878bandedResolutionMD5 = "2eacf3a41f71d7a035223512886dc2c7";
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.GVCF, WExIntervals, NA12878bandedResolutionMD5});
        tests.add(new Object[]{NA12878_WEx + " -I " + privateTestDir + "NA20313.highCoverageRegion.bam -sn NA12878",
                ReferenceConfidenceMode.GVCF, WExIntervals, NA12878bandedResolutionMD5});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.GVCF, WExIntervals + " -disableOptimizations", NA12878bandedResolutionMD5});

        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name = "MyDataProviderTetraploid")
    public Object[][] makeMyDataProviderTetraploid() {
        List<Object[]> tests = new ArrayList<>();

        final String PCRFreeIntervals = "-L 20:10,000,000-10,010,000";
        final String WExIntervals = "-L 20:10,000,000-10,100,000 -isr INTERSECTION -L " + hg19Chr20Intervals;

        // this functionality can be adapted to provide input data for whatever you might want in your data
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.NONE, PCRFreeIntervals, "c72403e4796942adfd272ecb28689df9"});
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.BP_RESOLUTION, PCRFreeIntervals, "3f446f0393a7d9eee75e107cafb1015b"});
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.GVCF, PCRFreeIntervals, "6cb4b227cb733cbcb888094c4465cf6a"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.NONE, WExIntervals, "e6efbcda94a9fda3d6a54e8ff7d22720"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.BP_RESOLUTION, WExIntervals, "01e7e4e268cb2d7ce50bac5dc30c7b62"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.GVCF, WExIntervals, "1d34a74bcec3513426f18e74f77a49e9"});

        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name = "MyDataProviderManyploid")
    public Object[][] makeMyDataProviderManyploid() {
        List<Object[]> tests = new ArrayList<>();

        final String PCRFreeIntervals = "-L 20:10,000,000-10,010,000";
        final String WExIntervals = "-L 20:10,000,000-10,100,000 -isr INTERSECTION -L " + hg19Chr20Intervals;

        // this functionality can be adapted to provide input data for whatever you might want in your data
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.NONE, PCRFreeIntervals, "47952bc713fb3c2d2beb8eca9ca4f2f2"});
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.BP_RESOLUTION, PCRFreeIntervals, "571931ace4b657d93ed295140963640b"});
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.GVCF, PCRFreeIntervals, "403c700ad4374b337802056e79373538"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.NONE, WExIntervals, "b699371662758d86d9d4b30f2b10f4c2"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.BP_RESOLUTION, WExIntervals, "7131959df133fafe2ac7884f7e09f687"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.GVCF, WExIntervals, "10c9a44f6933271dd64063607655abc2"});

        return tests.toArray(new Object[][]{});
    }

    /**
     * Test HaplotypeCaller, using MyDataProvider
     */
    @Test(dataProvider = "MyDataProvider")
    public void testHCWithGVCF(String bam, ReferenceConfidenceMode mode, String intervals, String md5) {
        final String commandLine = String.format("-T HaplotypeCaller --disableDithering --pcr_indel_model NONE -R %s -I %s %s -ERC %s --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, bam, intervals, mode, GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final String name = "testHCWithGVCF bam=" + bam + " intervals= " + intervals + " gvcf= " + mode;
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList(md5));
        executeTest(name, spec);
    }

    /**
     * Test HaplotypeCaller with haploid samples, using MyDataProviderHaploid
     */
    @Test(dataProvider = "MyDataProviderHaploid", enabled=true)
    public void testHCWithGVCFHaploid(final String bam, final ReferenceConfidenceMode mode, final String intervals, final String md5) {
        final String commandLine = String.format("-T HaplotypeCaller -ploidy 1 --disableDithering --pcr_indel_model NONE -R %s -I %s %s -ERC %s --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, bam, intervals, mode, GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final String name = "testHCWithGVCFHaploid bam=" + bam + " intervals= " + intervals + " gvcf= " + mode;
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList(md5));
        executeTest(name, spec);
    }

    /**
     * Test HaplotypeCaller with tetraploid samples, using MyDataProviderTetraploid
     */
    @Test(dataProvider = "MyDataProviderTetraploid", enabled=true)
    public void testHCWithGVCFTetraploid(final String bam, final ReferenceConfidenceMode mode, final String intervals, final String md5) {
        final String commandLine = String.format("-T HaplotypeCaller -ploidy 4 --disableDithering --pcr_indel_model NONE -R %s -I %s %s -ERC %s --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, bam, intervals, mode, GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final String name = "testHCWithGVCFTetraploid bam=" + bam + " intervals= " + intervals + " gvcf= " + mode;
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList(md5));
        executeTest(name, spec);
    }

    /**
     * Test HaplotypeCaller with manyploid samples, using MyDataProviderManyploid
     */
    @Test(dataProvider = "MyDataProviderManyploid", enabled=true)
    public void testHCWithGVCFManyploid(final String bam, final ReferenceConfidenceMode mode, final String intervals, final String md5) {
        final String commandLine = String.format("-T HaplotypeCaller -ploidy 33 --disableDithering --pcr_indel_model NONE -R %s -I %s %s -ERC %s --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, bam, intervals, mode, GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final String name = "testHCWithGVCFManyploid bam=" + bam + " intervals= " + intervals + " gvcf= " + mode;
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList(md5));
        executeTest(name, spec);
    }

    @Test
    public void testERCRegionWithNoCalledHaplotypes() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, privateTestDir + "noCallRefModel.bam", "20:17000001-18000001", GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList(""));
        spec.disableShadowBCF();
        executeTest("testERCRegionWithNoCalledHaplotypes", spec);
    }

    @Test()
    public void testMissingGVCFIndexException() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF",
                b37KGReference, privateTestDir + "noCallRefModel.bam", "20:17000001-18000001");
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", 1, UserException.GVCFIndexException.class);
        spec.disableShadowBCF();
        executeTest("testMissingGVCFIndexingStrategyException", spec);
    }

    /**
     * Test HaplotypeCaller to ensure it does not throw an exception when a .g.vcf output file is specified and the indexing arguments are omitted
     */
    @Test()
    public void testGVCFIndexNoThrow() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF",
                b37KGReference, privateTestDir + "noCallRefModel.bam", "20:17000000-17000100");
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine  + " -o %s", Arrays.asList(GATKVCFUtils.GVCF_EXT), Arrays.asList(""));
        spec.disableShadowBCF();
        executeTest("testGVCFIndexNoThrow", spec);
    }

    /**
     * Test HaplotypeCaller to ensure it does not throw an exception when a .g.vcf.gz output file is specified and the indexing arguments are omitted.
     * Verify that the output file is using the GZIP file format.
     */
    @Test()
    public void testGVCFGzIndexNoThrow() throws IOException {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF",
                b37KGReference, privateTestDir + "noCallRefModel.bam", "20:17000000-17000100");
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine, Arrays.asList(""));
        final File outputFile = createTempFile("testGVCFGzIndexNoThrow", "." + GATKVCFUtils.GVCF_GZ_EXT);
        spec.setOutputFileLocation(outputFile);
        spec.disableShadowBCF();
        executeTest("testGVCFGzIndexNoThrow", spec);
        final GZIPInputStream gzipOutputFileStream = new GZIPInputStream(new FileInputStream(outputFile));
    }

    @Test()
    public void testWrongParameterGVCFIndexException() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, privateTestDir + "noCallRefModel.bam", "20:17000001-18000001", GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER + 1);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", 1, UserException.GVCFIndexException.class);
        spec.disableShadowBCF();
        executeTest("testMissingGVCFIndexingStrategyException", spec);
    }

    @Test()
    public void testWrongTypeGVCFIndexException() {
        // ensure non-optimal, if optimal changes
        GATKVCFIndexType type = GATKVCFIndexType.DYNAMIC_SEEK;
        if (GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE == GATKVCFIndexType.DYNAMIC_SEEK)
            type = GATKVCFIndexType.DYNAMIC_SIZE;

        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, privateTestDir + "noCallRefModel.bam", "20:17000001-18000001", type, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", 1, UserException.GVCFIndexException.class);
        spec.disableShadowBCF();
        executeTest("testMissingGVCFIndexingStrategyException", spec);
    }

    private final static String WRONG_GVCF_RECORD_ORDER_BUGFIX_INTERVALS = privateTestDir + "gvcf_unsorted_records_bug.interval_list";
    private final static String WRONG_GVCF_RECORD_ORDER_BUGFIX_BAM = privateTestDir + "gvcf_unsorted_records_bug.bam";

    @Test()
    public void testWrongGVCFNonVariantRecordOrderBugFix() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, WRONG_GVCF_RECORD_ORDER_BUGFIX_BAM, WRONG_GVCF_RECORD_ORDER_BUGFIX_INTERVALS, GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList("5d13d8f0db7106e0873eaf1fc40db3ef"));
        spec.disableShadowBCF();
        executeTest("testMissingGVCFIndexingStrategyException", spec);
    }

    private static final String NOCALL_GVCF_BUGFIX_INTERVALS = privateTestDir + "gvcf_nocall_bug.interval_list";
    private static final String NOCALL_GVCF_BUGFIX_BAM = privateTestDir + "gvcf_nocall_bug.bam";
    private static final String GENERAL_PLOIDY_BUGFIX1_BAM = privateTestDir + "general-ploidy-arrayindex-bug-1.bam";
    private static final String GENERAL_PLOIDY_BUGFIX1_INTERVALS = privateTestDir + "general-ploidy-arrayindex-bug-1.intervals";
    private static final String GENERAL_PLOIDY_BUGFIX2_BAM = privateTestDir + "general-ploidy-arrayindex-bug-2.bam";
    private static final String GENERAL_PLOIDY_BUGFIX2_INTERVALS = privateTestDir + "general-ploidy-arrayindex-bug-2.intervals";


    @Test
    public void testNoCallGVCFMissingPLsBugFix() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, NOCALL_GVCF_BUGFIX_BAM, NOCALL_GVCF_BUGFIX_INTERVALS, GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList("c11f2be9a1dc4340b188b749c2b92929"));
        spec.disableShadowBCF();
        executeTest("testNoCallGVCFMissingPLsBugFix", spec);
    }

    /**
     * Checking that the bug's Exception is no longer been thrown (thus no md5).
     */
    @Test(enabled=true)
    public void testGeneralPloidyArrayIndexBug1Fix() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d -ploidy 1 -maxAltAlleles 2 -isr INTERSECTION -L 1:23696115-23696189",
                b37KGReference, GENERAL_PLOIDY_BUGFIX1_BAM, GENERAL_PLOIDY_BUGFIX1_INTERVALS, GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList(""));
        spec.disableShadowBCF();
        executeTest(" testGeneralPloidyArrayIndexBug1Fix", spec);
    }

    /**
     * Checking that the bug's Exception is no longer been thrown (thus no md5).
     */
    @Test(enabled=true)
    public void testGeneralPloidyArrayIndexBug2Fix() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d -ploidy 2 -maxAltAlleles 2 -A DepthPerSampleHC -A StrandBiasBySample -L 1:38052860-38052937",
                b37KGReference, GENERAL_PLOIDY_BUGFIX2_BAM, GENERAL_PLOIDY_BUGFIX2_INTERVALS, GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList(""));
        spec.disableShadowBCF();
        executeTest(" testGeneralPloidyArrayIndexBug2Fix", spec);
    }

    @Test
    public void testAlleleSpecificAnnotations() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d -G Standard -G AS_Standard --disableDithering",
                b37KGReference, privateTestDir + "NA12878.HiSeq.b37.chr20.10_11mb.bam", "20:10433000-10437000", GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Collections.singletonList("45347695ff437a9615163e9f3bcf8320"));
        spec.disableShadowBCF();
        executeTest(" testAlleleSpecificAnnotations", spec);
    }

    @Test
    public void testASMQMateRankSumAnnotation() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d -A AS_MQMateRankSumTest --disableDithering",
                b37KGReference, privateTestDir + "NA12878.HiSeq.b37.chr20.10_11mb.bam", "20:10433000-10437000", GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Collections.singletonList("5ca7ed6a7ff152b04bb2a7a97607091a"));
        spec.disableShadowBCF();
        executeTest(" testASMQMateRankSumAnnotation", spec);
    }

    @Test
    public void testBetaTestingAnnotationGroup() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d -G BetaTesting --disableDithering",
                b37KGReference, privateTestDir + "NA12878.HiSeq.b37.chr20.10_11mb.bam", "20:10433000-10437000", GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Collections.singletonList("a1677568cbd9304df927feea8a23085b"));
        spec.disableShadowBCF();
        executeTest(" testASMQMateRankSumAnnotation", spec);
    }

    @Test
    public void testASInsertSizeRankSum() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d -G Standard -G AS_Standard --disableDithering -A AS_InsertSizeRankSum",
                b37KGReference, privateTestDir + "NA12878.HiSeq.b37.chr20.10_11mb.bam", "20:10433000-10437000", GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Collections.singletonList("7809baf72434fdfa94e8e5862596e36b"));
        spec.disableShadowBCF();
        executeTest(" testASInsertSizeRankSum", spec);
    }

    @Test
    public void testHaplotypeCallerMultiAllelicNonRef() {
        final String commandLine = String.format("-T HaplotypeCaller -R %s -I %s -L %s -ERC GVCF --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d -A StrandAlleleCountsBySample",
                b37KGReference, privateTestDir + "multiallelic-nonref.bam", "2:47641259-47641859", GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Collections.singletonList("279cd69241f26eed2227ee22b568f333"));
        spec.disableShadowBCF();
        executeTest(" testHaplotypeCallerMultiAllelicNonRef", spec);
    }

    @Test
    public void testHaplotypeCallerMaxNumPLValues() {
        final String commandLine = String.format("-T HaplotypeCaller -R %s -I %s -L %s -ERC GVCF --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d -ploidy 4 -maxNumPLValues 70",
                b37KGReference, privateTestDir + "NA12878.HiSeq.b37.chr20.10_11mb.bam", validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.test.intervals", GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Collections.singletonList("29d2edbbc3c5dda8a9b502fa4337cf30"));
        spec.disableShadowBCF();
        executeTest("testHaplotypeCallerMaxNumPLValues", spec);
    }

    @Test
    public void testHaplotypeCallerMaxNumPLValuesExceededWithWarnLogLevel() throws IOException {
        // Need to see log WARN messages
        final Level level = logger.getLevel();
        logger.setLevel(Level.WARN);

        final File logFile = createTempFile("testMaxNumPLValuesExceededWithWarnLogLevel.log", ".tmp");
        final String logFileName = logFile.getAbsolutePath();

        final String commandLine = String.format("-T HaplotypeCaller -R %s -I %s -L %s -ERC GVCF --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d -ploidy 4 -maxNumPLValues 30 -log %s",
                b37KGReference, privateTestDir + "NA12878.HiSeq.b37.chr20.10_11mb.bam", validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.test.intervals",
                GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER, logFileName);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList("b2bf762651da36349319ece7628bb8aa"));
        spec.disableShadowBCF();
        executeTest("testHaplotypeCallerMaxNumPLValuesExceededWithWarnLogLevel", spec);
        // Make sure the "Maximum allowed number of PLs exceeded" messages are in the log
        Assert.assertTrue(FileUtils.readFileToString(logFile).contains("Maximum allowed number of PLs (30) exceeded for sample NA12878 at 20:10097101-10097101 with 35 possible genotypes. " +
                "No PLs will be output for these genotypes (which may cause incorrect results in subsequent analyses) unless the --max_num_PL_values argument is increased accordingly. " +
                "Unless the DEBUG logging level is used, this warning message is output just once per run and further warnings are suppressed."));
        Assert.assertFalse(FileUtils.readFileToString(logFile).contains("Maximum allowed number of PLs (30) exceeded for sample NA12878 at 20:10316239-10316241 with 70 possible genotypes."));
        Assert.assertTrue(FileUtils.readFileToString(logFile).contains("Maximum allowed number of PLs (30) was exceeded 2 time(s); the largest number of PLs found was 70."));
        // Set the log level back
        logger.setLevel(level);
    }

    @Test
    public void testHaplotypeCallerMaxNumPLValuesExceededWithDebugLogLevel() throws IOException {
        // Need to see log DEBUG messages
        final Level level = logger.getLevel();
        logger.setLevel(Level.DEBUG);

        final File logFile = createTempFile("testMaxNumPLValuesExceededWithDebugLogLevel.log", ".tmp");
        final String logFileName = logFile.getAbsolutePath();

        final String commandLine = String.format("-T HaplotypeCaller -R %s -I %s -L %s -ERC GVCF --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d -ploidy 4 -maxNumPLValues 30 -log %s",
                b37KGReference, privateTestDir + "NA12878.HiSeq.b37.chr20.10_11mb.bam", validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.test.intervals",
                GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER, logFileName);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList("b2bf762651da36349319ece7628bb8aa"));
        spec.disableShadowBCF();
        executeTest("testHaplotypeCallerMaxNumPLValuesExceededWithDebugLogLevel", spec);
        // Make sure the "Maximum allowed number of PLs exceeded" messages are in the log
        Assert.assertTrue(FileUtils.readFileToString(logFile).contains("Maximum allowed number of PLs (30) exceeded for sample NA12878 at 20:10097101-10097101 with 35 possible genotypes."));
        Assert.assertTrue(FileUtils.readFileToString(logFile).contains("Maximum allowed number of PLs (30) exceeded for sample NA12878 at 20:10316239-10316241 with 70 possible genotypes."));
        Assert.assertTrue(FileUtils.readFileToString(logFile).contains("Maximum allowed number of PLs (30) was exceeded 2 time(s); the largest number of PLs found was 70."));
        // Set the log level back
        logger.setLevel(level);
    }

    //Regression test for https://github.com/broadinstitute/gsa-unstable/issues/1345
    @Test
    public void testHaplotypeCallerGVCFBlocks() {
        final String commandLine = String.format("-T HaplotypeCaller -R %s -I %s -L 1:1-1000000 -ERC GVCF --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, privateTestDir + "gvcf_blocks_test.bam", GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList("a9ac40d37f4ac398162ce131a258a661"));
        spec.disableShadowBCF();
        executeTest("testHaplotypeCallerGVCFBlocks", spec);
    }

    @DataProvider(name = "dataBadGQBValues")
    public Object[][] dataBadGQBValues() {
        return new Object[][]{
                {Arrays.asList(-1, 10, 20)},
                {Arrays.asList(10, 20, 1)},
                {Arrays.asList(10, 10, 20)},
                {Arrays.asList(10, 20, VCFConstants.MAX_GENOTYPE_QUAL + 2)}
        };
    }
    @Test(dataProvider = "dataBadGQBValues")
    public void testBadGQBValues(final List<Integer> inputGQBValues) {
        final String inputGQBValuesString = inputGQBValues.stream().map(gqb -> "-GQB " + gqb).collect(Collectors.joining(" "));
        final String commandLine = String.format("-T HaplotypeCaller -R %s -I %s -L 1:1-1000000 -ERC GVCF %s --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, privateTestDir + "gvcf_blocks_test.bam", inputGQBValuesString, GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", 1, UserException.BadArgumentValue.class);
        spec.disableShadowBCF();
        executeTest("testBadGQBValues", spec);
    }

    @Test
    public void testHaplotypeCallerGVCSpanDel() {
        final String commandLine = String.format("-T HaplotypeCaller -R %s -I %s -L 1:26357667 -ERC GVCF --no_cmdline_in_header -A AS_ReadPosRankSumTest -A ReadPosRankSumTest -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, privateTestDir + "NexPond-377866-1:26357600-26357700.bam", GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList("5b1b75afcb1580df25a1cd7099ce050b"));
        spec.disableShadowBCF();
        executeTest("testHaplotypeCallerGVCSpanDel", spec);
    }
}
