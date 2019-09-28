package qut;

import java.util.concurrent.Callable;
import edu.au.jacobi.pattern.Match;

import static qut.Parallel.*;

public class DataBank implements Callable<Prediction> {

    public static int counterExec = 0;
    private final Gene gene;
    private final Gene referenceGene;
    private final NucleotideSequence nucleotides;

    DataBank(Gene gene, Gene refGene, NucleotideSequence nucleotides) {
        this.gene = gene;
        this.referenceGene = refGene;
        this.nucleotides = nucleotides;
    }

    @Override
    public Prediction call() throws Exception {
        if (Homologous(this.gene.sequence, this.referenceGene.sequence)) {
            NucleotideSequence upStreamRegion = GetUpstreamRegion(this.nucleotides, this.gene);
            Match prediction = PredictPromoter(upStreamRegion);
            if (prediction != null) {
                this.counterExec++;
                return new Prediction(this.referenceGene.name, prediction);
            }
        }
        return null;
    }

    public Gene getGene() {
        return this.gene;
    }

    public String getReferenceGeneName() {
        return this.referenceGene.name;
    }

    public PeptideSequence getReferenceGeneSeq() {
        return this.referenceGene.sequence;
    }

    public NucleotideSequence getNucleotides() {
        return this.nucleotides;
    }
}
