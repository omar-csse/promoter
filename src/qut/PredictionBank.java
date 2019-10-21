package qut;

import edu.au.jacobi.pattern.Match;

public class PredictionBank {

    private final String name;
    private final Match prediction;

    PredictionBank(String referenceGeneName, Match prediction) {
        this.name = referenceGeneName;
        this.prediction = prediction;
    }

    public String getName() {
        return this.name;
    }

    public Match getPrediction() {
        return this.prediction;
    }

}
