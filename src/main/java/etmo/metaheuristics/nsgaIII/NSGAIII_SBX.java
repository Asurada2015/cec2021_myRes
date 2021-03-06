//package etmo.metaheuristics.nsgaIII;
//
//import etmo.core.*;
//import etmo.util.JMException;
//import etmo.util.Ranking;
//import etmo.util.vector.TwoLevelWeightVectorGenerator;
//import etmo.util.vector.VectorGenerator;
//
//public class NSGAIII_SBX extends Algorithm {
//
//	private int populationSize_;
//
//	private int div1_;
//	private int div2_;
//
//	private SolutionSet population_;
//	SolutionSet offspringPopulation_;
//	SolutionSet union_;
//
//	int evaluations_;
//
//	Operator crossover_;
//	Operator mutation_;
//	Operator selection_;
//
//	double[][] lambda_; // reference points
//
//	boolean normalize_; // do normalization or not
//
//
//	public NSGAIII_SBX(ProblemSet problemSet) {
//		super(problemSet);
//	} // NSGAII
//
//	public SolutionSet execute() throws JMException, ClassNotFoundException {
//		int maxEvaluations_;
//
//		evaluations_ = 0;
//
//		maxEvaluations_ = ((Integer) this.getInputParameter("maxEvaluations"))
//				.intValue();
//
//		div1_ = ((Integer) this.getInputParameter("div1")).intValue();
//		div2_ = ((Integer) this.getInputParameter("div2")).intValue();
//
//
//		normalize_ = ((Boolean) this.getInputParameter("normalize")).booleanValue();
//
//		VectorGenerator vg = new TwoLevelWeightVectorGenerator(div1_, div2_,problemSet_.get(0).getNumberOfObjectives());
//		lambda_ = vg.getVectors();
//
//		populationSize_ = vg.getVectors().length;
//		if (populationSize_ % 2 != 0)
//			populationSize_ += 1;
//
//
//		mutation_ = operators_.get("mutation");
//		crossover_ = operators_.get("crossover");
//		selection_ = operators_.get("selection");
//
//		initPopulation();
//
//		while (evaluations_ < maxEvaluations_) {
//			offspringPopulation_ = new SolutionSet(populationSize_);
//			Solution[] parents = new Solution[2];
//			for (int i = 0; i < (populationSize_/2); i++) {
//				if (evaluations_ < maxEvaluations_) {
//					// obtain parents
//
//					parents = (Solution[]) selection_.execute(population_);
//
//					Solution[] offSpring = (Solution[]) crossover_.execute(parents);
//
//					mutation_.execute(offSpring[0]);
//					mutation_.execute(offSpring[1]);
//
//					problemSet_.get(0).evaluate(offSpring[0]);
//					problemSet_.get(0).evaluateConstraints(offSpring[0]);
//					problemSet_.get(0).evaluate(offSpring[1]);
//					problemSet_.get(0).evaluateConstraints(offSpring[1]);
//
//					offspringPopulation_.add(offSpring[0]);
//					offspringPopulation_.add(offSpring[1]);
//
//				} // if
//			} // for
//
//			union_ = ((SolutionSet) population_).union(offspringPopulation_);
//
//			// Ranking the union
//			Ranking ranking = new Ranking(union_);
//
//			int remain = populationSize_;
//			int index = 0;
//			SolutionSet front = null;
//			population_.clear();
//
//			// Obtain the next front
//			front = ranking.getSubfront(index);
//
//			while ((remain > 0) && (remain >= front.size())) {
//
//				for (int k = 0; k < front.size(); k++) {
//					population_.add(front.get(k));
//				} // for
//
//				// Decrement remain
//				remain = remain - front.size();
//
//				// Obtain the next front
//				index++;
//				if (remain > 0) {
//					front = ranking.getSubfront(index);
//				} // if
//			}
//
//			if (remain > 0) { // front contains individuals to insert
//
//				new Niching(population_, front, lambda_, remain, normalize_)
//						.execute();
//				remain = 0;
//			}
//
//			evaluations_=evaluations_+populationSize_;
//
//		}
//		Ranking ranking = new Ranking(population_);
//		return ranking.getSubfront(0);
//
//	}
//
//	public void initPopulation() throws JMException, ClassNotFoundException {
//
//		population_ = new SolutionSet(populationSize_);
//
//		for (int i = 0; i < populationSize_; i++) {
//			Solution newSolution = new Solution(problemSet_.get(0));
//
//			problemSet_.get(0).evaluate(newSolution);
//			problemSet_.get(0).evaluateConstraints(newSolution);
//
//			population_.add(newSolution);
//		} // for
//	} // initPopulation
//
//
//} // NSGA-III
//
