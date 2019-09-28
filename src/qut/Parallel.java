package qut;

import qut.*;
import jaligner.*;
import jaligner.matrix.*;
import edu.au.jacobi.pattern.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.Collectors;

public class Parallel {

    private static int numProcesses = Runtime.getRuntime().availableProcessors();
    private static HashMap<String, Sigma70Consensus> consensus = new HashMap<String, Sigma70Consensus>();
    private static ThreadLocal<Series> sigma70_pattern = ThreadLocal.withInitial(() -> Sigma70Definition.getSeriesAll_Unanchored(0.7));
    private static final Matrix BLOSUM_62 = BLOSUM62.Load();
    private static byte[] complement = new byte['z'];

    static {
        complement['C'] = 'G'; complement['c'] = 'g';
        complement['G'] = 'C'; complement['g'] = 'c';
        complement['T'] = 'A'; complement['t'] = 'a';
        complement['A'] = 'T'; complement['a'] = 't';
    }


    private static List<Gene> ParseReferenceGenes(String referenceFile) throws FileNotFoundException, IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(referenceFile)));
        List<Gene> referenceGenes = new ArrayList<Gene>();
        while (true) {
            String name = reader.readLine();
            if (name == null)
                break;
            String sequence = reader.readLine();
            referenceGenes.add(new Gene(name, 0, 0, sequence));
            consensus.put(name, new Sigma70Consensus());
        }
        consensus.put("all", new Sigma70Consensus());
        reader.close();
        return referenceGenes;
    }

    public static boolean Homologous(PeptideSequence A, PeptideSequence B) {
        return SmithWatermanGotoh.align(new Sequence(A.toString()), new Sequence(B.toString()), BLOSUM_62, 10f, 0.5f).calculateScore() >= 60;
    }

    public static NucleotideSequence GetUpstreamRegion(NucleotideSequence dna, Gene gene) {
        int upStreamDistance = 250;
        if (gene.location < upStreamDistance)
            upStreamDistance = gene.location-1;

        if (gene.strand == 1)
            return new NucleotideSequence(java.util.Arrays.copyOfRange(dna.bytes, gene.location-upStreamDistance-1, gene.location-1));
        else
        {
            byte[] result = new byte[upStreamDistance];
            int reverseStart = dna.bytes.length - gene.location + upStreamDistance;
            for (int i=0; i<upStreamDistance; i++)
                result[i] = complement[dna.bytes[reverseStart-i]];
            return new NucleotideSequence(result);
        }
    }

    public static Match PredictPromoter(NucleotideSequence upStreamRegion) {
        return BioPatterns.getBestMatch(sigma70_pattern.get(), upStreamRegion.toString());
    }

    private static void ProcessDir(List<String> list, File dir) {
        if (dir.exists()) {
            for (File file : dir.listFiles()) {
                if (file.isDirectory()) {
                    ProcessDir(list, file);
                } else {
                    if ( !file.getPath().contains(".DS_Store") )
                        list.add(file.getPath());
                }
            }
        }
    }

    private static List<String> ListGenbankFiles(String dir) throws IOException {
        List<String> list = new ArrayList<String>();
        ProcessDir(list, new File(dir));
        List<String> result = list.stream().filter(path -> !path.contains(".DS_Store")).collect(Collectors.toList());
        return result;
    }

    private static GenbankRecord Parse(String file) throws IOException {
        GenbankRecord record = new GenbankRecord();
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
        record.Parse(reader);
        reader.close();
        return record;
    }

    public static void run(String referenceFile, String dir) throws FileNotFoundException, IOException, InterruptedException, ExecutionException {

        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        List<DataBank> dataBanks = new LinkedList<>();

        for (String filename : ListGenbankFiles(dir)) {
            System.out.println(filename);
            GenbankRecord record = Parse(filename);
            for (Gene referenceGene : referenceGenes) {
                System.out.println(referenceGene.name);
                for (Gene gene : record.genes) {
                    dataBanks.add(new DataBank(gene, referenceGene, record.nucleotides));
                }
            }
        }

        ExecutorService executer = Executors.newFixedThreadPool(numProcesses);
        List<Future<Prediction>> results = executer.invokeAll(dataBanks);

        for (Future<Prediction> future: results) {
            Prediction prediction = future.get();
            if (prediction != null) {
                consensus.get(prediction.getName()).addMatch(prediction.getPrediction());
                consensus.get("all").addMatch(prediction.getPrediction());
            }
        }

        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
            System.out.println(entry.getKey() + " " + entry.getValue());

        executer.shutdown();
    }

    public static void main(String[] args) throws FileNotFoundException, IOException, InterruptedException, ExecutionException {
        String currentdir = System.getProperty("user.dir");
        double startTime = System.currentTimeMillis();
        run(currentdir+"/referenceGenes.list", currentdir+"/Ecoli");
        double endTime = System.currentTimeMillis();
        System.out.println("\n\nOverall Time (m): " + ( (endTime - startTime) / 1000 / 60 ) + "m");
        System.out.println("Overall Time (s): " + ( (endTime - startTime) / 1000 ) + "s");
    }
}
