package se.kth.warrenk;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.commons.cli.*;
import org.jgrapht.Graph;
import org.slf4j.Logger;
import uk.ac.ox.well.cortexjdk.Main;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
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

/**
 * Created by kiran on 11/09/2017.
 */
public class NovelExonFinder {
    public static void main(String[] args) throws Exception {
        // Parse command line
        Options o = new Options();
        o.addOption("g",  "graph",                true, "The joined transcriptome graph");
        o.addOption("s",  "seeds",                true, "The FASTA file representing seed sequences");
        o.addOption("tn", "transcriptomeName",    true, "The sample name assigned to the canonical transcriptome in the graph");
        o.addOption("tl", "transcriptomeLinksDb", true, "The linksdb file for the canonical transcriptome");
        o.addOption("sn", "sampleName",           true, "The sample name assigned to the sample of interest in the graph");
        o.addOption("sl", "sampleLinksDb",        true, "The linksdb file for the sample");
        o.addOption("o",  "output",               true, "The output file");

        CommandLineParser clp = new DefaultParser();
        CommandLine cl = clp.parse(o, args);

        if (!cl.hasOption("graph") ||
            !cl.hasOption("transcriptomeLinksDb") ||
            !cl.hasOption("sampleLinksDb") ||
            !cl.hasOption("seeds") ||
            !cl.hasOption("output") ||
            !cl.hasOption("transcriptomeName") ||
            !cl.hasOption("sampleName")) {

            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("java -jar novelexonfinder.jar", o);

            System.exit(1);
        }

        CortexGraph cg = new CortexGraph(cl.getOptionValue("graph"));
        CortexLinks rl = new CortexLinks(cl.getOptionValue("transcriptomeLinksDb"));
        CortexLinks al = new CortexLinks(cl.getOptionValue("sampleLinksDb"));
        FastaSequenceFile mads = new FastaSequenceFile(new File(cl.getOptionValue("seeds")), true);
        PrintStream out = new PrintStream(new File(cl.getOptionValue("output")));

        int transcriptomeColor = cg.getColorForSampleName(cl.getOptionValue("transcriptomeName"));
        int sampleColor = cg.getColorForSampleName(cl.getOptionValue("sampleName"));

        Logger log = Main.getLogger();

        // Loads mads sequences
        Set<String> madsKmers = new HashSet<>();
        ReferenceSequence rseq;
        int numMadsSeqs = 0;
        while ((rseq = mads.nextSequence()) != null) {
            String seq = rseq.getBaseString();
            for (int i = 0; i <= rseq.length() - cg.getKmerSize(); i++) {
                String sk = seq.substring(i, i + cg.getKmerSize());
                madsKmers.add(sk);
            }

            numMadsSeqs++;
        }
        log.info("Loaded {} mads sequences into {} kmers", numMadsSeqs, madsKmers.size());

        // Look for novel kmers
        Set<CortexKmer> novelKmers = new HashSet<>();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing kmers")
                .message("kmers processed")
                .maxRecord(cg.getNumRecords())
                .make(log);

        for (CortexRecord cr : cg) {
            if (isNovel(cr, transcriptomeColor, sampleColor)) {
                novelKmers.add(cr.getCortexKmer());
            }

            pm.update();
        }

        log.info("  found {} novel kmers", novelKmers.size());

        // Traverse contigs at novel kmers
        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(sampleColor)
                .graph(cg)
                .links(rl, al)
                .make();

        pm = new ProgressMeterFactory()
                .header("Examining mads kmers")
                .message("novel mads examined")
                .maxRecord(madsKmers.size())
                .make(log);

        Set<CortexKmer> processedKmers = new TreeSet<>();
        for (String madsKmer : madsKmers) {
            CortexKmer ck = new CortexKmer(madsKmer);

            if (!processedKmers.contains(ck)) {
                List<CortexVertex> w = e.walk(madsKmer);

                boolean found = false;
                for (CortexVertex v : w) {
                    if (!found && novelKmers.contains(v.getCk())) {
                        found = true;
                    }

                    processedKmers.add(v.getCk());
                }

                if (found) {
                    log.info("  found seed={} contig_length={} contig={}", madsKmer, w.size(), TraversalEngine.toContig(w));

                    out.println(">" + madsKmer);
                    out.println(TraversalEngine.toContig(w));
                }
            }

            pm.update();
        }
    }

    private static boolean isNovel(CortexRecord cr, int transcriptomeColor, int sampleColor) {
        boolean transcriptomeHasCoverage = cr.getCoverage(transcriptomeColor) > 0;
        boolean sampleHasCoverage = cr.getCoverage(sampleColor) > 0;

        return !transcriptomeHasCoverage && sampleHasCoverage;
    }
}
