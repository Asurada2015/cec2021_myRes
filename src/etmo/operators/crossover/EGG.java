//  SBXCrossover.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package etmo.operators.crossover;

import jmetal.core.Solution;
import jmetal.encodings.solutionType.ArrayRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.operators.crossover.Crossover;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.wrapper.XReal;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

/**
 * This class allows to apply a SBX crossover operator using two parent
 * solutions.
 */
public class EGG extends Crossover {
  /**
   * EPS defines the minimum difference allowed between real values
   */
  private static final double EPS= 1.0e-14;
                                                                                      
  private static final double ETA_C_DEFAULT_ = 20.0;
  private Double crossoverProbability_ = 0.9 ;
  private double distributionIndex_ = ETA_C_DEFAULT_;
  private double F_ = 0.5;
  
  private double Dr_ = 0.2;
  private double Lr_ = 0.2;
  private double Er_ = 0.7;

  /**
   * Valid solution types to apply this operator 
   */
  private static final List VALID_TYPES = Arrays.asList(RealSolutionType.class,
                                                  ArrayRealSolutionType.class) ;
  
  /** 
   * Constructor
   * Create a new SBX crossover operator whit a default
   * index given by <code>DEFAULT_INDEX_CROSSOVER</code>
   */
  public EGG(HashMap<String, Object> parameters) {
  	super (parameters) ;
  	
  	if (parameters.get("probability") != null)
  		crossoverProbability_ = (Double) parameters.get("probability") ;  		
  	if (parameters.get("distributionIndex") != null)
  		distributionIndex_    = (Double) parameters.get("distributionIndex") ; 
  	if (parameters.get("Dr") != null)
  		Dr_ = (Double) parameters.get("Dr");  
  	if (parameters.get("Lr") != null)
  		Lr_ = (Double) parameters.get("Lr");  
  	if (parameters.get("Er") != null)
  		Er_ = (Double) parameters.get("Er");
  } // SBXCrossover
    
  /**
   * Perform the crossover operation. 
   * @param probability Crossover probability
   * @param parent1 The first parent
   * @param parent2 The second parent
   * @return An array containing the two offsprings
   */
  public Solution doCrossover(double probability,
                              Solution parent1,
                              Solution parent2,
                              Solution parent3) throws JMException {
    
    Solution offSpring = new Solution();

	offSpring = new Solution(parent1);
                    
	int i;
    double rand;
    double yL, yu;
    double c1, c2, c3;
    double alpha, beta, betaq;
    double valueX1,valueX2,valueX3;
	XReal x1 = new XReal(parent1);
	XReal x2 = new XReal(parent2);
	XReal x3 = new XReal(parent3);
	XReal offs = new XReal(offSpring);
	int numberOfVariables = x1.getNumberOfDecisionVariables() ;
	for (i=0; i<numberOfVariables; i++){
		valueX1 = x1.getValue(i);
		valueX2 = x2.getValue(i);
        valueX3 = x3.getValue(i);
        yL = x1.getLowerBound(i) ;
        yu = x1.getUpperBound(i) ;
        
//        if (PseudoRandom.randDouble()<= P1_ ){//StdRandom.gaussian(0.3, 0.02) Dr=Disburbance rate
//        	if (java.lang.Math.abs(valueX1- valueX2) > EPS){
        		rand = PseudoRandom.randDouble();
	            alpha = 2.0;
	            if (rand <= (1.0/alpha)){
	              betaq = Math.pow ((rand*alpha),(1.0/(distributionIndex_+1.0)));
	            } else {
	              betaq = Math.pow ((1.0/(2.0 - rand*alpha)),(1.0/(distributionIndex_+1.0)));
	            } // if
	            c1 = 0.5*((1+betaq)*valueX1+(1-betaq)*valueX2);
	            rand = PseudoRandom.randDouble();
	            alpha = 2.0;
	            if (rand <= (1.0/alpha)){
	              betaq = Math.pow ((rand*alpha),(1.0/(distributionIndex_+1.0)));
	            } else {
	              betaq = Math.pow ((1.0/(2.0 - rand*alpha)),(1.0/(distributionIndex_+1.0)));
	            } // if
	            c2 = 0.5*((1-betaq)*valueX1+(1+betaq)*valueX2);
	            c3 = valueX1 + F_ * (valueX2-valueX3);
            
	            if (c1<yL) c1=yL;
	            if (c2<yL) c2=yL;
	            if (c1>yu) c1=yu;
	            if (c2>yu) c2=yu;
	            if (c3<yL) c3=yL;
	            if (c3>yu) c3=yu;
	            double r= PseudoRandom.randDouble();
	            if(r<=Dr_) {
		            if (PseudoRandom.randDouble() <= Lr_ ){
		            	offs.setValue(i,c3);
		            }else{
		            	if (PseudoRandom.randDouble() <= 0.5 ){
				            offs.setValue(i, c2);
				        } else {
				            offs.setValue(i, c1);
				        }
		            }// if
	            } else if (r <= Dr_+Er_){ //
		        	offs.setValue(i, valueX2);
	        	}else{
	        		offs.setValue(i, valueX1);
	        	}
	} // for                       
	return offSpring;                                                                                      
  } // doCrossover
  
  
  /**
  * Executes the operation
  * @param object An object containing an array of two parents
  * @return An object containing the offSprings
  */
  public Object execute(Object object) throws JMException {
    Solution[] parents = (Solution[])object;
    	
    if (parents.length != 3) {
      Configuration.logger_.severe("SBXCrossover.execute: operator needs two " +
          "parents");
      Class cls = String.class;
      String name = cls.getName(); 
      throw new JMException("Exception in " + name + ".execute()") ;
    } // if

    if (!(VALID_TYPES.contains(parents[0].getType().getClass())  &&
          VALID_TYPES.contains(parents[1].getType().getClass())  &&
          VALID_TYPES.contains(parents[2].getType().getClass())) ) {
      Configuration.logger_.severe("SBXCrossover.execute: the solutions " +
					"type " + parents[0].getType() + " is not allowed with this operator");

      Class cls = String.class;
      String name = cls.getName(); 
      throw new JMException("Exception in " + name + ".execute()") ;
    } // if 
    
    if (parameters_.get("probability") != null)
  		crossoverProbability_ = (Double) parameters_.get("probability") ;  		
  	if (parameters_.get("distributionIndex") != null)
  		distributionIndex_    = (Double) parameters_.get("distributionIndex") ;
  	if (parameters_.get("F") != null)
  		F_ = (Double) parameters_.get("F") ;  	
  	if (parameters_.get("Dr") != null)
  		Dr_ = (Double) parameters_.get("Dr") ;  
  	if (parameters_.get("Lr") != null)
  		Lr_ = (Double) parameters_.get("Lr") ;  
  	if (parameters_.get("Er") != null)
  		Er_ = (Double) parameters_.get("Er") ;
    Solution offSpring;
    offSpring = doCrossover(crossoverProbability_,
                            parents[0],
                            parents[1],
                            parents[2]);
    return offSpring;
  } // execute 
} // SBXCrossover
