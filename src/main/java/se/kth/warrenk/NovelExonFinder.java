package se.kth.warrenk;

import org.jgrapht.Graph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.io.File;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

/**
 * Created by kiran on 11/09/2017.
 */
public class NovelExonFinder {
    public static void main(String[] args) throws Exception {
        if (args.length != 4) { throw new RuntimeException("Wrong number of args"); }

        String graphFile = args[0];
        String refLinks = args[1];
        String altLinks = args[2];
        String outFile = args[3];

        int transcriptomeColor = 0;
        int madsColor = 2;

        CortexGraph cg = new CortexGraph(graphFile);
        CortexLinks rl = new CortexLinks(refLinks);
        CortexLinks al = new CortexLinks(altLinks);

        PrintStream out = new PrintStream(new File(outFile));

        Set<CortexRecord> novelKmers = new HashSet<>();

        System.err.println("Processing kmers...");
        for (CortexRecord cr : cg) {
            if (isNovel(cr, transcriptomeColor, madsColor)) {
                novelKmers.add(cr);
            }
        }
        System.err.println("\t" + novelKmers.size() + " novel kmers found");

        TraversalEngine e = new TraversalEngineFactory()
                .joiningColors(madsColor)
                .graph(cg)
                .traversalDirection(BOTH)
                .combinationOperator(OR)
                .stoppingRule(JoiningStopper.class)
                .links(rl, al)
                .make();

        Set<CortexKmer> processedKmers = new TreeSet<>();
        for (CortexRecord nk : novelKmers) {
            if (!processedKmers.contains(nk.getCortexKmer())) {
                for (int c = 0; c < nk.getNumColors(); c++) {
                    if (c != transcriptomeColor && c != madsColor && nk.getCoverage(c) > 0) {
                        e.getConfiguration().setTraversalColor(c);

                        List<CortexVertex> w = e.walk(nk.getKmerAsString());

                        boolean touchedMadsGene = false;
                        for (CortexVertex cv : w) {
                            if (cv.getCr().getCoverage(madsColor) > 0) {
                                touchedMadsGene = true;
                            }

                            processedKmers.add(cv.getCk());
                        }

                        if (touchedMadsGene) {
                            String contig = TraversalEngine.toContig(w);

                            out.println(">" + nk.getKmerAsString());
                            out.println(contig);
                        }
                    }
                }
            }
        }
    }

    private static boolean isNovel(CortexRecord cr, int transcriptomeColor, int madsColor) {
        boolean transcriptomeHasCoverage = cr.getCoverage(transcriptomeColor) > 0;
        boolean madsboxHasCoverage = cr.getCoverage(madsColor) > 0;

        if (!transcriptomeHasCoverage && !madsboxHasCoverage) {
            for (int c = 0; c < cr.getNumColors(); c++) {
                if (c != transcriptomeColor && c != madsColor && cr.getCoverage(c) > 0) {
                    return true;
                }
            }
        }

        return false;
    }
}
