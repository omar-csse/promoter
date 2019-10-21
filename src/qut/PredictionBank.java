package qut;

import edu.au.jacobi.pattern.Match;

public class Prediction {

    private final String name;
    private final Match prediction;

    Prediction(String referenceGeneName, Match prediction) {
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
