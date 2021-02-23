package etmo.metaheuristics.maoeac;

import etmo.core.Algorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.SolutionSet;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.problems.benchmarks_ETMO.*;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.Ranking;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;

public class MaOEAC_main {
    public static void main(String[] args) throws IOException, JMException, ClassNotFoundException {
        ProblemSet problemSet1; // The problem to solve
        ProblemSet problemSet2;
        Algorithm algorithm; // The algorithm to use
        Operator crossover; // Crossover operator
        Operator mutation; // Mutation operator
        Operator selection;

        HashMap parameters; // Operator parameters

        for (int pCase = 9; pCase <= 16; pCase++ ){
            switch (pCase){
                case 1:
                    problemSet1 = ETMOF1.getProblem();
                    break;
                case 2:
                    problemSet1 = ETMOF2.getProblem();
                    break;
                case 3:
                    problemSet1 = ETMOF3.getProblem();
                    break;
                case 4:
                    problemSet1 = ETMOF4.getProblem();
                    break;
                case 5:
                    problemSet1 = ETMOF5.getProblem();
                    break;
                case 6:
                    problemSet1 = ETMOF6.getProblem();
                    break;
                case 7:
                    problemSet1 = ETMOF7.getProblem();
                    break;
                case 8:
                    problemSet1 = ETMOF8.getProblem();
                    break;
                case 9:
                    problemSet1 = ETMOF9.getProblem();
                    break;
                case 10:
                    problemSet1 = ETMOF10.getProblem();
                    break;
                case 11:
                    problemSet1 = ETMOF11.getProblem();
                    break;
                case 12:
                    problemSet1 = ETMOF12.getProblem();
                    break;
                case 13:
                    problemSet1 = ETMOF13.getProblem();
                    break;
                case 14:
                    problemSet1 = ETMOF14.getProblem();
                    break;
                case 15:
                    problemSet1 = ETMOF15.getProblem();
                    break;
                case 16:
                    problemSet1 = ETMOF16.getProblem();
                    break;
                default:
                    problemSet1 = ETMOF1.getProblem();
            }


            int taskNumber = problemSet1.size();
//            System.out.println("taskNumber = "+taskNumber);
            System.out.println("TaskID\t" + "IGD for " + problemSet1.get(0).getName()+" to " +problemSet1.get(taskNumber-1).getName());
            for (int tsk=0;tsk<taskNumber;tsk++) {

                problemSet2 = problemSet1.getTask(tsk);
                algorithm = new MaOEAC(problemSet2);

                String pf = "PF/StaticPF/" + problemSet2.get(0).getHType() + "_" + problemSet2.get(0).getNumberOfObjectives() + "D.pf";
                //String pf = "PF/StaticPF/" + "convex.pf";
                //System.out.println(pf);
                algorithm.setInputParameter("populationSize", 100);
                algorithm.setInputParameter("maxGenerations",1000);

                parameters = new HashMap();
                parameters.put("probability", 1.0);
                parameters.put("distributionIndex", 20.0);
//                crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);
                crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

                // Mutation operator
                parameters = new HashMap();
                parameters.put("probability", 1.0 / problemSet2.getMaxDimension());
                parameters.put("distributionIndex", 20.0);
                mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

                // Selection Operator
                parameters = null ;
                selection = SelectionFactory.getSelectionOperator("BinaryTournament", parameters) ;

                // Add the operators to the algorithm
                algorithm.addOperator("crossover", crossover);
                algorithm.addOperator("mutation", mutation);
                algorithm.addOperator("selection", selection);

//                System.out.println("RunID\t" + "IGD for " + problemSet2.get(0).getName());
                DecimalFormat form = new DecimalFormat("#.####E0");
                QualityIndicator indicator = new QualityIndicator(problemSet2.get(0), pf);
                int times = 21;
                double aveIGD = 0;
                for (int i = 1; i <= times; i++) {
                    SolutionSet population = algorithm.execute();
                    Ranking ranking = new Ranking(population);
                    population = ranking.getSubfront(0);
//                    population.printObjectivesToFile("MaOEAC_"+problemSet2.get(0).getNumberOfObjectives()+"Obj_"+
//                            problemSet2.get(0).getName()+ "_" + problemSet2.get(0).getNumberOfVariables() + "D_run"+i+".txt");
                    double igd = indicator.getIGD(population);
                    aveIGD += igd;
//                    System.out.println(i + "\t" + form.format(igd));
                }
//                System.out.println("Average IGD for " + problemSet2.get(0).getName() + ": " + form.format(aveIGD / times));
//                System.out.println();
//                System.out.println(tsk+"\t"+form.format(aveIGD / times));
                System.out.println(form.format(aveIGD / times));

            }
        }



    }
}
