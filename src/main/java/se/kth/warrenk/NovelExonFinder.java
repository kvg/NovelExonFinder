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
        if (args.length != 3) { throw new RuntimeException("Wrong number of args"); }

        String graphFile = args[0];
        String refLinks = args[1];
        String altLinks = args[2];

        int transcriptomeColor = 0;
        int madsColor = 2;

        CortexGraph cg = new CortexGraph(graphFile);
        CortexLinks rl = new CortexLinks(refLinks);
        CortexLinks al = new CortexLinks(altLinks);

        Set<CortexRecord> novelKmers = new HashSet<>();

        System.out.println("Processing kmers...");
        for (CortexRecord cr : cg) {
            if (isNovel(cr, transcriptomeColor, madsColor)) {
                novelKmers.add(cr);
            }
        }
        System.out.println("\t" + novelKmers.size() + " novel kmers found");

        TraversalEngine e = new TraversalEngineFactory()
                .joiningColors(madsColor)
                .graph(cg)
                .traversalDirection(BOTH)
                .combinationOperator(OR)
                .stoppingRule(JoiningStopper.class)
                .links(rl, al)
                .make();

        Set<CortexKmer> seen = new TreeSet<>();
        for (CortexRecord nk : novelKmers) {
            if (!seen.contains(nk.getCortexKmer())) {
                for (int c = 0; c < nk.getNumColors(); c++) {
                    if (c != transcriptomeColor && c != madsColor && nk.getCoverage(c) > 0) {
                        e.getConfiguration().setTraversalColor(c);

                        //List<CortexVertex> w = e.walk(nk.getKmerAsString());
                        Graph<CortexVertex, CortexEdge> g = e.dfs(nk.getKmerAsString());

                        if (g != null && g.vertexSet().size() > 0) {
                            System.out.println(">" + nk.getKmerAsString());

                            for (CortexVertex cv : g.vertexSet()) {
                                seen.add(cv.getCk());
                            }
                        }


                        /*
                        String contig = TraversalEngine.toContig(w);


                        System.out.println(">" + nk.getKmerAsString());
                        System.out.println(contig);
                        */
                    }
                }
            }
        }
    }

    private static boolean isNovel(CortexRecord cr, int transcriptomeColor, int madsColor) {
        boolean transcriptomeHasCoverage = cr.getCoverage(transcriptomeColor) > 0;

        if (!transcriptomeHasCoverage) {
            for (int c = 0; c < cr.getNumColors(); c++) {
                if (c != transcriptomeColor && c != madsColor && cr.getCoverage(c) > 0) {
                    return true;
                }
            }
        }

        return false;
    }
}
