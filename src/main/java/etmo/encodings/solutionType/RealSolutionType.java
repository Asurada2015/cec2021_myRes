//  RealSolutionType.java
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

package etmo.encodings.solutionType;

import etmo.core.ProblemSet;
import etmo.core.SolutionType;
import etmo.core.Variable;
import etmo.encodings.variable.Real;
import etmo.util.JMException;

/**
 * Class representing a solution type composed of real variables
 */
public class RealSolutionType extends SolutionType {

	/**
	 * Constructor
	 * 
	 * @param problemSet
	 *            Problem to solve
	 */
	public RealSolutionType(ProblemSet problemSet) {
		super(problemSet);
	} // Constructor

	/**
	 * Creates the variables of the solution
	 */
	public Variable[] createVariables() throws JMException {
		Variable[] variables = new Variable[problemSet_.getMaxDimension()];

		for (int var = 0; var < problemSet_.getMaxDimension(); var++)
			variables[var] = new Real(problemSet_.getUnifiedLowerLimit(), problemSet_.getUnifiedUpperLimit());

//		for (int var = 0; var < problemSet_.getMaxDimension(); var++){
//			variables[var] = new Real(problemSet_.getUnifiedLowerLimit(), problemSet_.getUnifiedUpperLimit());
//			variables[var].setValue(0.4);
//		}




		return variables;
	} // createVariables
} // RealSolutionType
