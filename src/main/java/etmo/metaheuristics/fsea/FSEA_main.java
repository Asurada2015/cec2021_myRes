package etmo.metaheuristics.fsea;

import etmo.core.Algorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.SolutionSet;
import etmo.metaheuristics.nsgaII.NSGAII;
import etmo.metaheuristics.utils.printIGD;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.problems.benchmarks_CEC2017.*;
import etmo.problems.benchmarks_ETMO.*;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;

public class FSEA_main {
    public static void main(String args[]) throws IOException, JMException, ClassNotFoundException {
        ProblemSet problemSet1; // The problem to solve
        ProblemSet problemSet2;
        Algorithm algorithm; // The algorithm to use
        Operator crossover; // Crossover operator
        Operator mutation; // Mutation operator
        Operator selection;

        HashMap parameters; // Operator parameters
        for (int pCase = 1; pCase <= 8; pCase++ ){
            switch (pCase){
//                case 1:
//                    problemSet1 = CIHS.getProblem();
//                    break;
//                case 2:
//                    problemSet1 = CIMS.getProblem();
//                    break;
//                case 3:
//                    problemSet1 = CILS.getProblem();
//                    break;
//                case 4:
//                    problemSet1 = PIHS.getProblem();
//                    break;
//                case 5:
//                    problemSet1 = PIMS.getProblem();
//                    break;
//                case 6:
//                    problemSet1 = PILS.getProblem();
//                    break;
//                case 7:
//                    problemSet1 = NIHS.getProblem();
//                    break;
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
            int times = 21;
            double cpIGD[][] = new double[taskNumber][times];

            for (int tsk = 0; tsk < taskNumber; tsk++) {
//            for (int tsk=0; tsk < taskNumber; tsk++) {

                problemSet2 = problemSet1.getTask(tsk);
                algorithm = new FSEA(problemSet2);

//                String pf = "PF/StaticPF/" + problemSet2.get(0).getHType() + "_" + problemSet2.get(0).getNumberOfObjectives() + "D.csv";
                String pfcal = "PF/StaticPF/" + problemSet2.get(0).getHType() + "_" + problemSet2.get(0).getNumberOfObjectives() + "D.pf";


//                add chart
//                List<Double> referencePoint = new ArrayList<>() ;



//                ChartContainerWithReferencePoints chart = new ChartContainerWithReferencePoints("test", 80);
//                chart.setFrontChart(0,1, pf);
//                chart.setReferencePoint(convertReferencePointListToListOfLists(referencePoint, problemSet2.getNumberOfObjs(0)));
//                chart.initChart();

                algorithm.setInputParameter("populationSize", 100);
                algorithm.setInputParameter("maxEvaluations", 100 * 1000);

                parameters = new HashMap();
                parameters.put("probability", 0.9);
                parameters.put("distributionIndex", 20.0);
                crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

                // Mutation operator
                parameters = new HashMap();
                parameters.put("probability", 1.0 / problemSet2.getMaxDimension());
                parameters.put("distributionIndex", 20.0);
                mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

                // Selection Operator
                parameters = null ;
                selection = SelectionFactory.getSelectionOperator("BinaryTournament2", parameters) ;

                // Add the operators to the algorithm
                algorithm.addOperator("crossover", crossover);
                algorithm.addOperator("mutation", mutation);
                algorithm.addOperator("selection", selection);

//                System.out.println("RunID\t" + "IGD for " + problemSet2.get(0).getName());
                DecimalFormat form = new DecimalFormat("#.####E0");
                QualityIndicator indicator = new QualityIndicator(problemSet2.get(0), pfcal);
//				int times = 21;
                double aveIGD = 0;
                for (int i = 1; i <= times; i++) {
                    SolutionSet population = algorithm.execute();
//                    chart.saveChart("moead", BitmapEncoder.BitmapFormat.PNG);
//                    population.printObjectivesToFile("MOEAD_"+problemSet2.get(0).getNumberOfObjectives()+"Obj_"+
//                        problemSet2.get(0).getName()+ "_" + problemSet2.get(0).getNumberOfVariables() + "D_run"+i+".txt");
                    double igd = indicator.getIGD(population);
                    aveIGD += igd;
                    cpIGD[tsk][i - 1] = igd;

//                    System.out.println(i + "\t" + form.format(igd));
//                    System.out.println(i + "\t" + form.format(igd));
//                    resources/weightVectorFiles/moead/W3D_100.dat


                }
//                System.out.println("Average IGD for " + problemSet2.get(0).getName() + ": " + form.format(aveIGD / times));
                System.out.println(form.format(aveIGD / times));
                //                System.out.println();
            }
//            String path = "NSGAII_2017F1-7.txt";
//            printIGD.printIGDtoText(path, cpIGD, taskNumber, times);


        }

//		follow the origin


//		problemSet1 = ETMOF5.getProblem();
//		int taskNumber = problemSet1.size();
//		System.out.println("taskNumber = "+taskNumber);
//		for (int tsk=0;tsk<taskNumber;tsk++) {
//
//			problemSet2 = problemSet1.getTask(tsk);
//			algorithm = new NSGAII(problemSet2);
//
//			String pf = "PF/StaticPF/" + problemSet2.get(0).getHType() + "_" + problemSet2.get(0).getNumberOfObjectives() + "D.pf";
//			//String pf = "PF/StaticPF/" + "convex.pf";
//		    //System.out.println(pf);
//			algorithm.setInputParameter("populationSize", 100);
//			algorithm.setInputParameter("maxEvaluations", 100 * 1000);
//
//			parameters = new HashMap();
//			parameters.put("probability", 0.9);
//			parameters.put("distributionIndex", 20.0);
//			crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);
//
//			// Mutation operator
//			parameters = new HashMap();
//			parameters.put("probability", 1.0 / problemSet2.getMaxDimension());
//			parameters.put("distributionIndex", 20.0);
//			mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);
//
//			// Selection Operator
//		    parameters = null ;
//		    selection = SelectionFactory.getSelectionOperator("BinaryTournament2", parameters) ;
//
//			// Add the operators to the algorithm
//			algorithm.addOperator("crossover", crossover);
//			algorithm.addOperator("mutation", mutation);
//			algorithm.addOperator("selection", selection);
//
//			System.out.println("RunID\t" + "IGD for " + problemSet2.get(0).getName());
//			DecimalFormat form = new DecimalFormat("#.####E0");
//			QualityIndicator indicator = new QualityIndicator(problemSet2.get(0), pf);
//			int times = 1;
//			double aveIGD = 0;
//				for (int i = 1; i <= times; i++) {
//					SolutionSet population = algorithm.execute();
//					Ranking ranking = new Ranking(population);
//					population = ranking.getSubfront(0);
//					population.printObjectivesToFile("NSGAII_"+problemSet2.get(0).getNumberOfObjectives()+"Obj_"+
//					                                 problemSet2.get(0).getName()+ "_" + problemSet2.get(0).getNumberOfVariables() + "D_run"+i+".txt");
//					double igd = indicator.getIGD(population);
//					aveIGD += igd;
//					System.out.println(i + "\t" + form.format(igd));
//				}
//				System.out.println("Average IGD for " + problemSet2.get(0).getName() + ": " + form.format(aveIGD / times));
//				System.out.println();
//			}
    }
}
