package se.kth.warrenk;

import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;

/**
 * Created by kiran on 11/09/2017.
 */
public class NovelExonFinder {
    public static void main(String[] args) throws Exception {
        if (args.length < 2) { throw new RuntimeException("Not enough args"); }

        String graphFile = args[0];
        int transcriptomeColor = Integer.valueOf(args[1]);

        CortexGraph cg = new CortexGraph(graphFile);

        TraversalEngine e = new TraversalEngineFactory()
                .joiningColors(transcriptomeColor)
                .graph(cg)
                .make();

        Set<CortexRecord> novelKmers = new HashSet<>();

        System.out.println("Processing kmers...");
        for (CortexRecord cr : cg) {
            if (isNovel(cr, transcriptomeColor)) {
                novelKmers.add(cr);
            }
        }
        System.out.println("\t" + novelKmers.size() + " novel kmers found");
    }

    private static boolean isNovel(CortexRecord cr, int transcriptomeColor) {
        boolean transcriptomeHasCoverage = cr.getCoverage(transcriptomeColor) > 0;

        if (!transcriptomeHasCoverage) {
            for (int c = 0; c < cr.getNumColors(); c++) {
                if (c != transcriptomeColor && cr.getCoverage(c) > 0) {
                    return true;
                }
            }
        }

        return false;
    }
}
