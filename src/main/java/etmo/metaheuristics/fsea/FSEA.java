package etmo.metaheuristics.fsea;

import etmo.core.*;
import etmo.util.Distance;
import etmo.util.JMException;
import etmo.util.Ranking;
import etmo.util.comparators.CrowdingComparator;

public class FSEA extends Algorithm {
    /**
     * Constructor
     *
     * @param problemSet The problem to be solved
     */
    public FSEA(ProblemSet problemSet) {
        super(problemSet);
    }

    @Override
    public SolutionSet execute() throws JMException, ClassNotFoundException {
        int populationSize;
        int maxEvaluations;
        int evaluations;


        SolutionSet population;
        SolutionSet offspringPopulation;
        SolutionSet union;

        Operator mutationOperator;
        Operator crossoverOperator;
        Operator selectionOperator;

        Distance distance = new Distance();

        // Read the parameters
        populationSize = ((Integer) getInputParameter("populationSize")).intValue();
        maxEvaluations = ((Integer) getInputParameter("maxEvaluations")).intValue();

        // Initialize the variables
        population = new SolutionSet(populationSize);
        evaluations = 0;


        // Read the operators
        mutationOperator = operators_.get("mutation");
        crossoverOperator = operators_.get("crossover");
        selectionOperator = operators_.get("selection");

        // Create the initial solutionSet
        Solution newSolution;
        for (int i = 0; i < populationSize; i++) {
            newSolution = new Solution(problemSet_);
            problemSet_.get(0).evaluate(newSolution);
            problemSet_.get(0).evaluateConstraints(newSolution);
            evaluations++;
            population.add(newSolution);
        } // for

        // Generations
        while (evaluations < maxEvaluations) {

            // Create the offSpring solutionSet
            offspringPopulation = new SolutionSet(populationSize);
            Solution[] parents = new Solution[2];
            for (int i = 0; i < (populationSize / 2); i++) {
                if (evaluations < maxEvaluations) {
                    // obtain parents
                    parents[0] = (Solution) selectionOperator.execute(population);
                    parents[1] = (Solution) selectionOperator.execute(population);
                    Solution[] offSpring = (Solution[]) crossoverOperator.execute(parents);
                    mutationOperator.execute(offSpring[0]);
                    mutationOperator.execute(offSpring[1]);
                    problemSet_.get(0).evaluate(offSpring[0]);
                    problemSet_.get(0).evaluateConstraints(offSpring[0]);
                    problemSet_.get(0).evaluate(offSpring[1]);
                    problemSet_.get(0).evaluateConstraints(offSpring[1]);
                    offspringPopulation.add(offSpring[0]);
                    offspringPopulation.add(offSpring[1]);
                    evaluations += 2;
                } // if
            } // for

            // Create the solutionSet union of solutionSet and offSpring
            union = ((SolutionSet) population).union(offspringPopulation);

            // Ranking the union
            Ranking ranking = new Ranking(union);


            int remain = populationSize;
            int index = 0;
            SolutionSet front = null;
            population.clear();

            // Obtain the next front
            front = ranking.getSubfront(index);
//            System.out.println(evaluations + ": ranking nums is " + front.size());
//            for (int i = 0; i < front.size(); i++){
//                Solution a = front.get(i);
//                Variable[] v = a.getDecisionVariables();
////                for (int j = 0; j < v.length; j++){
////                    System.out.printf("%f ", v[j].getValue());
////                }
//                System.out.println(v[0].getValue());
//            }





            while ((remain > 0) && (remain >= front.size())) {
                // Assign crowding distance to individuals
                distance.crowdingDistanceAssignment(front, problemSet_.get(0).getNumberOfObjectives());
                // Add the individuals of this front
                for (int k = 0; k < front.size(); k++) {
                    population.add(front.get(k));
                } // for

                // Decrement remain
                remain = remain - front.size();

                // Obtain the next front
                index++;
                if (remain > 0) {
                    front = ranking.getSubfront(index);
                } // if
            } // while

            // Remain is less than front(index).size, insert only the best one
            if (remain > 0) { // front contains individuals to insert
                distance.crowdingDistanceAssignment(front, problemSet_.get(0).getNumberOfObjectives());
                front.sort(new CrowdingComparator());
                for (int k = 0; k < remain; k++) {
                    population.add(front.get(k));
                } // for

                remain = 0;
            } // if

        } // while

        Ranking ranking = new Ranking(population);

        return ranking.getSubfront(0);
    }
}
