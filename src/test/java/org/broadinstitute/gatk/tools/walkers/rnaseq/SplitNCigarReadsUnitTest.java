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

package org.broadinstitute.gatk.tools.walkers.rnaseq;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.clipping.ReadClipperTestUtils;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

/**
 *
 * Tests all possible (and valid) cigar strings that might contain any cigar elements. It uses a code that were written to test the ReadClipper walker.
 * For valid cigar sting in length 8 there are few thousands options, with N in every possible option and with more than one N (for example 1M1N1M1N1M1N2M).
 * The cigarElements array is used to provide all the possible cigar element that might be included.
 *
 * User: ami
 * Date: 11/14/13
 * Time: 6:49 PM
 */
public class SplitNCigarReadsUnitTest extends BaseTest {
    final static CigarElement[] cigarElements = {
            new CigarElement(1, CigarOperator.HARD_CLIP),
            new CigarElement(1, CigarOperator.SOFT_CLIP),
            new CigarElement(1, CigarOperator.INSERTION),
            new CigarElement(1, CigarOperator.DELETION),
            new CigarElement(1, CigarOperator.MATCH_OR_MISMATCH),
            new CigarElement(1, CigarOperator.SKIPPED_REGION)
    };

    private CachingIndexedFastaSequenceFile referenceReader;
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void setup() throws FileNotFoundException {
        referenceReader = new CachingIndexedFastaSequenceFile(new File(hg18Reference));
        genomeLocParser = new GenomeLocParser(referenceReader.getSequenceDictionary());
    }

    private final class TestManager extends OverhangFixingManager {
        public TestManager() {
            super(null, genomeLocParser, referenceReader, 10000, 1, 40, false);
        }
    }

    @Test(enabled = true)
    public void splitReadAtN() {
        final int cigarStringLength = 10;
        final List<Cigar> cigarList = ReadClipperTestUtils.generateCigarList(cigarStringLength,cigarElements);

        // For Debugging use those lines (instead of above cigarList) to create specific read:
        //------------------------------------------------------------------------------------
        // final GATKSAMRecord tmpRead = GATKSAMRecord.createRandomRead(6);
        // tmpRead.setCigarString("1M1N1M");

        // final List<Cigar> cigarList = new ArrayList<>();
        // cigarList.add(tmpRead.getCigar());

        for(Cigar cigar: cigarList){

            final int numOfSplits = numOfNElements(cigar.getCigarElements());

            if(numOfSplits != 0 && isCigarDoesNotHaveEmptyRegionsBetweenNs(cigar)){

                final TestManager manager = new TestManager();
                GATKSAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar, 0);
                SplitNCigarReads.splitNCigarRead(read, manager);
                List<OverhangFixingManager.SplitRead> splitReads = manager.getReadsInQueueForTesting();
                final int expectedReads = numOfSplits+1;
                Assert.assertEquals(splitReads.size(),expectedReads,"wrong number of reads after split read with cigar: "+cigar+" at Ns [expected]: "+expectedReads+" [actual value]: "+splitReads.size());
                final List<Integer> readLengths = consecutiveNonNElements(read.getCigar().getCigarElements());
                int index = 0;
                int offsetFromStart = 0;
                for(final OverhangFixingManager.SplitRead splitRead: splitReads){
                    int expectedLength = readLengths.get(index);
                    Assert.assertTrue(splitRead.read.getReadLength() == expectedLength,
                            "the "+index+" (starting with 0) split read has a wrong length.\n" +
                            "cigar of original read: "+cigar+"\n"+
                            "expected length: "+expectedLength+"\n"+
                            "actual length: "+splitRead.read.getReadLength()+"\n");
                    assertBases(splitRead.read.getReadBases(), read.getReadBases(), offsetFromStart);
                    index++;
                    offsetFromStart += expectedLength;
                }
            }
        }
    }

    private int numOfNElements(final List<CigarElement> cigarElements){
        int numOfNElements = 0;
        for (CigarElement element: cigarElements){
            if (element.getOperator() == CigarOperator.SKIPPED_REGION)
                numOfNElements++;
        }
        return numOfNElements;
    }

    private static boolean isCigarDoesNotHaveEmptyRegionsBetweenNs(final Cigar cigar) {
        boolean sawM = false;
        boolean sawS = false;

        for (CigarElement cigarElement : cigar.getCigarElements()) {
            if (cigarElement.getOperator().equals(CigarOperator.SKIPPED_REGION)) {
                if(!sawM && !sawS)
                    return false;
                sawM = false;
                sawS = false;
            }
            if (cigarElement.getOperator().equals(CigarOperator.MATCH_OR_MISMATCH))
                sawM = true;
            if (cigarElement.getOperator().equals(CigarOperator.SOFT_CLIP))
                sawS = true;

        }
        if(!sawS && !sawM)
            return false;
        return true;
    }

    private List<Integer> consecutiveNonNElements(final List<CigarElement> cigarElements){
        final LinkedList<Integer> results = new LinkedList<>();
        int consecutiveLength = 0;
        for(CigarElement element: cigarElements){
            final CigarOperator op = element.getOperator();
            if(op.equals(CigarOperator.MATCH_OR_MISMATCH) || op.equals(CigarOperator.SOFT_CLIP) || op.equals(CigarOperator.INSERTION)){
                consecutiveLength += element.getLength();
            }
            else if(op.equals(CigarOperator.SKIPPED_REGION))
            {
                if(consecutiveLength != 0){
                    results.addLast(consecutiveLength);
                    consecutiveLength = 0;
                }
            }
        }
        if(consecutiveLength != 0)
            results.addLast(consecutiveLength);
        return results;
    }

    private void assertBases(final byte[] actualBase, final byte[] expectedBase, final int startIndex) {
        for (int i = 0; i < actualBase.length; i++) {
            Assert.assertEquals(actualBase[i], expectedBase[startIndex + i],"unmatched bases between: "+ Arrays.toString(actualBase)+"\nand:\n"+Arrays.toString(expectedBase)+"\nat position: "+i);
        }
    }

}
