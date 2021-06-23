package etmo.metaheuristics.moead;

import etmo.core.Algorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.SolutionSet;
import etmo.metaheuristics.utils.Plot2D_demo1;
import etmo.metaheuristics.utils.printIGD;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.problems.benchmarks_CEC2017.*;
import etmo.problems.benchmarks_CEC2019.*;
import etmo.problems.benchmarks_ETMO.*;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import org.knowm.xchart.BitmapEncoder;
import org.uma.jmetal.util.chartcontainer.ChartContainerWithReferencePoints;
import org.uma.jmetal.util.front.impl.ArrayFront;


import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class

MOEAD_main {
    public static void main(String args[]) throws IOException, JMException, ClassNotFoundException, InterruptedException {
        ProblemSet problemSet1; // The problem to solve
        ProblemSet problemSet2;

        Algorithm algorithm; // The algorithm to use
        Operator crossover; // Crossover operator
        Operator mutation; // Mutation operator
        Operator selection;

        HashMap parameters; // Operator parameters
        for (int pCase = 1; pCase <=10; pCase++ ){
            switch (pCase){
                case 1:
                    problemSet1 = CPLX1.getProblem();
                    break;
                case 2:
                    problemSet1 = CPLX2.getProblem();
                    break;
                case 3:
                    problemSet1 = CPLX3.getProblem();
                    break;
                case 4:
                    problemSet1 = CPLX4.getProblem();
                    break;
                case 5:
                    problemSet1 = CPLX5.getProblem();
                    break;
                case 6:
                    problemSet1 = CPLX6.getProblem();
                    break;
                case 7:
                    problemSet1 = CPLX7.getProblem();
                    break;
                case 8:
                    problemSet1 = CPLX8.getProblem();
                    break;
                case 9:
                    problemSet1 = CPLX9.getProblem();
                    break;
                case 10:
                    problemSet1 = CPLX10.getProblem();
                    break;

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
//                case 8:
//                    problemSet1 = NIMS.getProblem();
//                    break;
//                case 9:
//                    problemSet1 = NILS.getProblem();
//                    break;


//                case 1:
//                    problemSet1 = ETMOF1.getProblem();
//                    break;
//                case 2:
//                    problemSet1 = ETMOF2.getProblem();
//                    break;
//                case 3:
//                    problemSet1 = ETMOF4.getProblem();
//                    break;
//                case 4:
//                    problemSet1 = ETMOF5.getProblem();
//                    break;
//                case 5:
//                    problemSet1 = ETMOF6.getProblem();
//                    break;
//                case 6:
//                    problemSet1 = ETMOF7.getProblem();
//                    break;
//                case 7:
//                    problemSet1 = ETMOF8.getProblem();
//                    break;
//                case 8:
//                    problemSet1 = ETMOF3.getProblem();
//                    break;
//                case 9:
//                    problemSet1 = ETMOF9.getProblem();
//                    break;
//                case 10:
//                    problemSet1 = ETMOF10.getProblem();
//                    break;
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
                case 17:
                    problemSet1 = ETMOF17.getProblem();
                    break;
                case 18:
                    problemSet1 = ETMOF18.getProblem();
                    break;
                case 19:
                    problemSet1 = ETMOF19.getProblem();
                    break;
                case 20:
                    problemSet1 = ETMOF20.getProblem();
                    break;
                case 21:
                    problemSet1 = ETMOF21.getProblem();
                    break;
                case 22:
                    problemSet1 = ETMOF22.getProblem();
                    break;
                case 23:
                    problemSet1 = ETMOF23.getProblem();
                    break;
                case 24:
                    problemSet1 = ETMOF24.getProblem();
                    break;
                default:
                    problemSet1 = ETMOF1.getProblem();
            }

            int taskNumber = problemSet1.size();
            int times = 21;
//            System.out.println("taskNumber = "+taskNumber);
            double cpIGD[][] = new double[taskNumber][times];

            for (int tsk = 0; tsk < taskNumber; tsk++) {
//            for (int tsk=0; tsk < taskNumber; tsk++) {

                problemSet2 = problemSet1.getTask(tsk);
                algorithm = new MOEAD(problemSet2);

//                String pf = "PF/StaticPF/" + problemSet2.get(0).getHType() + "_" + problemSet2.get(0).getNumberOfObjectives() + "D.csv";
//                String pfcal = "PF/StaticPF/" + problemSet2.get(0).getHType() + "_" + problemSet2.get(0).getNumberOfObjectives() + "D.pf";
                String pfcal =  "PF/cec2019/" + problemSet2.get(0).getHType() + ".pf";
//
//
//                add chart
//                List<Double> referencePoint = new ArrayList<>() ;



//                ChartContainerWithReferencePoints chart = new ChartContainerWithReferencePoints("test", 80);
//                chart.setFrontChart(0,1, pf);
//                chart.setReferencePoint(convertReferencePointListToListOfLists(referencePoint, problemSet2.getNumberOfObjs(0)));
//                chart.initChart();

                algorithm.setInputParameter("populationSize", 10000);
                algorithm.setInputParameter("maxEvaluations", 100 * 100);

                algorithm.setInputParameter("dataDirectory", "resources/weightVectorFiles/moead");


                algorithm.setInputParameter("T", 20);
                algorithm.setInputParameter("delta", 0.9);
                algorithm.setInputParameter("nr", 2);

                parameters = new HashMap();
                parameters.put("CR", 1.0);
                parameters.put("F", 0.9);
                crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover",parameters);

                // Mutation operator
                parameters = new HashMap();
                parameters.put("probability", 1.0 / problemSet2.get(0).getNumberOfVariables());
                parameters.put("distributionIndex", 20.0);
                mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);


                algorithm.addOperator("crossover", crossover);
                algorithm.addOperator("mutation", mutation);

//                System.out.println("RunID\t" + "IGD for " + problemSet2.get(0).getName());
                DecimalFormat form = new DecimalFormat("#.####E0");
                QualityIndicator indicator = new QualityIndicator(problemSet2.get(0), pfcal);
//                int times = 21;
                double aveIGD = 0;
//                double cpIGD[][] = new double[taskNumber][times];
                for (int i = 1; i <= times; i++) {
                    SolutionSet population = algorithm.execute();
//                    chart.saveChart("moead", BitmapEncoder.BitmapFormat.PNG);
//                    population.printObjectivesToFile("MOEAD_"+problemSet2.get(0).getNumberOfObjectives()+"Obj_"+
//                        problemSet2.get(0).getName()+ "_" + problemSet2.get(0).getNumberOfVariables() + "D_run"+i+".txt");
//                    Plot2D_demo1 test = new Plot2D_demo1(pfcal);
//                    test.doPlot2D(population);

                    double igd = indicator.getIGD(population);
                    aveIGD += igd;
//                    System.out.println(i + "\t" + form.format(igd));
//                    System.out.println(i + "\t" + form.format(igd));
//                    resources/weightVectorFiles/moead/W3D_100.dat
                    cpIGD[tsk][i - 1] = igd;
                }
//                System.out.println("Average IGD for " + problemSet2.get(0).getName() + ": " + form.format(aveIGD / times));
                System.out.println(form.format(aveIGD / times));
                //                System.out.println();
            }
//            String path = "MOEAD_B_CEC2017.txt";
//            printIGD.printIGDtoText(path, cpIGD, taskNumber, times);

        }





    }

    private static List<List<Double>> convertReferencePointListToListOfLists(List<Double> referencePoints, int numberOfObjectives) {
        List<List<Double>> referencePointList;
        referencePointList = new ArrayList<>();

        for (int i = 0; i <= (referencePoints.size() - numberOfObjectives); i+=numberOfObjectives) {
            List<Double> newReferencePoint = new ArrayList<>(numberOfObjectives) ;
            for (int j = i; j < (i + numberOfObjectives); j++) {
                newReferencePoint.add(referencePoints.get(j)) ;
            }

            referencePointList.add(newReferencePoint) ;
        }

        return referencePointList ;
    }

}
