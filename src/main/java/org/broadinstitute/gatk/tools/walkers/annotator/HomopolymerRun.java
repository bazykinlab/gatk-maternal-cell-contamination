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

package org.broadinstitute.gatk.tools.walkers.annotator;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.GenomeLoc;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Largest contiguous homopolymer run of the variant allele
 *
 * <p>Repetitive sequences such as homopolymers are difficult to map to the reference because they are associated with multiple alignment possibilities. The proximity of a long homopolymer to your variant site increases the chance that reads were mapped incorrectly in the surrounding region and lowers confidence in the call. If there is a homopolymer on either side of a site, this annotation outputs the length of its largest run.</p>
 *
 * <h3>Caveats</h3>
 * <ul>
 *     <li>This can only be computed for bi-allelic sites.</li>
 *     <li>The calculation only looks at direct runs of the alternate allele adjacent to this position, which is not a very accurate method.</li>
 *     <li>This is an experimental annotation. As such, it is unsupported; we do not make any guarantees that it will work properly, and you use it at your own risk.</li>
 * </ul>
 */
public class HomopolymerRun extends InfoFieldAnnotation implements ExperimentalAnnotation {

    private final static Logger logger = Logger.getLogger(HomopolymerRun.class);

    private boolean ANNOTATE_INDELS = true;
    
    private final static int BAD_RUN_LENGTH = -1;

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {

        if ( !vc.isBiallelic() )
            return null;

        int run;
        if ( vc.isSNP() ) {
            run = computeHomopolymerRun(vc.getAlternateAllele(0).getBases()[0], ref);
        } else if ( vc.isIndel() && ANNOTATE_INDELS ) {
            run = computeIndelHomopolymerRun(vc,ref);
        } else {
            return null;
        }

        if ( run == BAD_RUN_LENGTH )
            return null;
        
        Map<String, Object> map = new HashMap<>();
        map.put(getKeyNames().get(0), String.format("%d", run));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList(GATKVCFConstants.HOMOPOLYMER_RUN_KEY); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0))); }

    public boolean useZeroQualityReads() { return false; }

    private static int computeHomopolymerRun(byte altAllele, ReferenceContext ref) {

        // TODO -- this needs to be computed in a more accurate manner
        // We currently look only at direct runs of the alternate allele adjacent to this position

        byte[] bases = ref.getBases();
        GenomeLoc window = ref.getWindow();
        GenomeLoc locus = ref.getLocus();

        int refBasePos = (int)(locus.getStart() - window.getStart());

        int leftRun = 0;
        for ( int i = refBasePos - 1; i >= 0; i--) {
            if ( bases[i] != altAllele )
                break;
            leftRun++;
        }

        int rightRun = 0;
        for ( int i = refBasePos + 1; i < bases.length; i++) {
            if ( bases[i] != altAllele )
                break;
            rightRun++;
        }

        return Math.max(leftRun, rightRun);
     }

    private static int computeIndelHomopolymerRun(VariantContext vc, ReferenceContext ref) {
        byte[] bases = ref.getBases();
        GenomeLoc locus = ref.getLocus();
        GenomeLoc window = ref.getWindow();
        int refBasePos = (int) (locus.getStart() - window.getStart())+1;
        if ( vc.isSimpleDeletion() ) {
            if ( refBasePos + vc.getReference().length() - 1 >= bases.length ) {
                logger.warn("Encountered a homopolymer at " + locus.toString() + " longer than the tool's default window size, so the position was skipped. " +
                        "To process this position, add --reference_window_stop to your command with a value equal or greater than " + vc.getReference().length());
                return BAD_RUN_LENGTH;
            }

            return computeHomopolymerRun(bases[refBasePos], ref);
        } else {
            // check that inserted bases are the same
            byte insBase = vc.getAlternateAllele(0).getBases()[0];
            for ( byte b : vc.getAlternateAllele(0).getBases() ) {
                if ( insBase != (char) b ) {
                    return 0;
                }
            }

            return computeHomopolymerRun(insBase,ref);
        }
    }
}