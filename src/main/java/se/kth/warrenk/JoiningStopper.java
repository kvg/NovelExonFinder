package se.kth.warrenk;

import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.AbstractTraversalStoppingRule;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

import java.util.Set;

/**
 * Created by kiran on 12/09/2017.
 */
public class JoiningStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    private int numJoinedKmersSeen = 0;

    @Override
    public boolean hasTraversalSucceeded(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph, DeBruijnGraph rois) {
        boolean hasJoined = false;

        for (int c : joiningColors) {
            hasJoined |= cv.getCr().getCoverage(c) > 0;
        }

        if (hasJoined) {
            numJoinedKmersSeen++;
        }

        return numJoinedKmersSeen >= 5;
    }

    @Override
    public boolean hasTraversalFailed(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph, DeBruijnGraph rois) {
        return numAdjacentEdges != 1;
    }
}
