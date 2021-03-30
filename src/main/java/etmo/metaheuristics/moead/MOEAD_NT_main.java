package etmo.metaheuristics.moead;

import etmo.core.*;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.problems.benchmarks_ETMO.*;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;

public class MOEAD_NT_main {
    public static void main(String[] args) throws IOException, JMException, ClassNotFoundException {
        ProblemSet problemSet;

        Algorithm algorithm; // The algorithm to use
        Operator crossover; // Crossover operator
        Operator mutation; // Mutation operator

        HashMap parameters; // Operator parameters
        for (int pCase = 1; pCase <= 7; pCase++ ){
            switch (pCase){
                case 1:
                    problemSet = ETMOF1.getProblem();
                    break;
                case 2:
                    problemSet = ETMOF2.getProblem();
                    break;
                case 3:
                    problemSet = ETMOF4.getProblem();
                    break;
                case 4:
                    problemSet = ETMOF5.getProblem();
                    break;
                case 5:
                    problemSet = ETMOF6.getProblem();
                    break;
                case 6:
                    problemSet = ETMOF7.getProblem();
                    break;
                case 7:
                    problemSet = ETMOF8.getProblem();
                    break;
                case 8:
                    problemSet = ETMOF8.getProblem();
                    break;
                case 9:
                    problemSet = ETMOF9.getProblem();
                    break;
                case 10:
                    problemSet = ETMOF10.getProblem();
                    break;
                case 11:
                    problemSet = ETMOF11.getProblem();
                    break;
                case 12:
                    problemSet = ETMOF12.getProblem();
                    break;
                case 13:
                    problemSet = ETMOF13.getProblem();
                    break;
                case 14:
                    problemSet = ETMOF14.getProblem();
                    break;
                case 15:
                    problemSet = ETMOF15.getProblem();
                    break;
                case 16:
                    problemSet = ETMOF16.getProblem();
                    break;
                default:
                    problemSet = ETMOF1.getProblem();
            }

            int taskNumber = problemSet.size();

            String[] pf = new String[taskNumber];
            for (int i = 0; i < pf.length; i++){
                pf[i] = "PF/StaticPF/" + problemSet.get(i).getHType() + "_" + problemSet.get(i).getNumberOfObjectives() + "D.pf";
            }

            algorithm = new MOEAD_NT(problemSet);

            algorithm.setInputParameter("populationSize", 100);
            algorithm.setInputParameter("maxEvaluations", 100 * taskNumber * 1000);

            algorithm.setInputParameter("dataDirectory", "resources/weightVectorFiles/moead");


            algorithm.setInputParameter("T", 20);
            algorithm.setInputParameter("delta", 0.9);
            algorithm.setInputParameter("nr", 2);

            parameters = new HashMap();
            parameters.put("CR", 1.0);
            parameters.put("F", 0.5);
            crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover",parameters);

            // Mutation operator
            parameters = new HashMap();
            parameters.put("probability", 1.0 / problemSet.get(0).getNumberOfVariables());
            parameters.put("distributionIndex", 20.0);
            mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);


            algorithm.addOperator("crossover", crossover);
            algorithm.addOperator("mutation", mutation);



//                System.out.println("RunID\t" + "IGD for " + problemSet2.get(0).getName());
            DecimalFormat form = new DecimalFormat("#.####E0");


            int times = 21;
            double ave[] = new double[taskNumber];

            for (int t = 1; t <= times; t++){
                SolutionSet population = algorithm.execute();

                SolutionSet[] resPopulation = new SolutionSet[taskNumber];
                for (int i = 0; i < problemSet.size(); i++)
                    resPopulation[i] = new SolutionSet();

                for (int i = 0; i < population.size(); i++) {
                    Solution sol = population.get(i);

                    int pid = i / (population.size() / taskNumber);

                    int start = problemSet.get(pid).getStartObjPos();
                    int end = problemSet.get(pid).getEndObjPos();

                    Solution newSolution = new Solution(end - start + 1);

                    for (int k = start; k <= end; k++)
                        newSolution.setObjective(k - start, sol.getObjective(k));

                    resPopulation[pid].add(newSolution);
                }

                double igd;
//				System.out.print(t + "\t");
                for(int i = 0; i < taskNumber; i++){
                    QualityIndicator indicator = new QualityIndicator(problemSet.get(i), pf[i]);
                    if(resPopulation[i].size()==0)
                        continue;
//				getTask中用到add影响problem起始和结束值
//				resPopulation[i].printObjectivesToFile("MOMFEA_"+problemSet.getTask(i).get(0).getNumberOfObjectives()+"Obj_"+
//						problemSet.getTask(i).get(0).getName()+ "_" + problemSet.getTask(i).get(0).getNumberOfVariables() + "D_run"+t+".txt");
//					resPopulation[i].printObjectivesToFile("MOMFEA_"+problemSet.get(i).getNumberOfObjectives()+"Obj_"+
//							problemSet.get(i).getName()+ "_" + problemSet.get(i).getNumberOfVariables() + "D_run"+t+".txt");

                    igd =  indicator.getIGD(resPopulation[i]);
//					System.out.print(form.format(igd) + "\t" );
                    ave[i] += igd;
                }



            }

            for(int i=0;i<taskNumber;i++)
//				System.out.println("Average IGD for " + problemSet.get(i).getName()+ ": " + form.format(ave[i] / times));
                System.out.println(form.format(ave[i] / times));

//            QualityIndicator indicator = new QualityIndicator(problemSet.get(0), pf[0]);
//            int times = 21;
//            double aveIGD = 0;
//            for (int i = 1; i <= times; i++) {
//                SolutionSet population = algorithm.execute();
////                    chart.saveChart("moead", BitmapEncoder.BitmapFormat.PNG);
////                    population.printObjectivesToFile("MOEAD_"+problemSet2.get(0).getNumberOfObjectives()+"Obj_"+
////                        problemSet2.get(0).getName()+ "_" + problemSet2.get(0).getNumberOfVariables() + "D_run"+i+".txt");
//                double igd = indicator.getIGD(population);
//                aveIGD += igd;
////                    System.out.println(i + "\t" + form.format(igd));
////                    System.out.println(i + "\t" + form.format(igd));
////                    resources/weightVectorFiles/moead/W3D_100.dat
//
//
//            }
////                System.out.println("Average IGD for " + problemSet2.get(0).getName() + ": " + form.format(aveIGD / times));
//            System.out.println(form.format(aveIGD / times));
//            //                System.out.println();


        }
    }
}
