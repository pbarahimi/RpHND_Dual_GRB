import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;

public class Main {

	private static double alpha = 0.2;
	private static double[][] tmpFlows = MyArray.read("w.txt");
	private static double[][] coordinates = MyArray.read("coordinates.txt");
	private static double[][] distances = Distance.get(coordinates);
	private static int nVar = tmpFlows.length;
	private static double[][] flows = new double[nVar][nVar];
	private static int D = 1; // Maximum number of simultaneous disruptions
	private static int R = (int) Math.pow(2, D + 1) - 2; // Largest index in the
															// full binary tree.
	private static int _R = (int) Math.pow(2, D) - 2; // Total number of nodes
														// in the solution tree excluding leaf-nodes.
	private static double q = 0.05;
	private static int P = 4; // number of hubs to be located
	private static int M = nVar * R; // the big M

	/**
	 * 
	 * @param i
	 * @param j
	 * @param k
	 * @param m
	 * @return operating probability of a route
	 */
	public static double Q(int i, int k, int m, int j) {
		double result = q;
		if (k != i && j != m)
			result = q + q(k, m);
		else if (m != j)
			result = q(i, m);
		else if (k != i)
			result = q(j, k);
		else if (i == k && j == m)
			result = 0;
		else
			System.out.println("Not include in the Q(i,k,m,j)!");
		return result;
	}

	/**
	 * Cikmj
	 * 
	 */
	private static double Cikmj(int i, int k, int m, int j) {
		double cost = distances[i][k] + (1 - alpha) * distances[k][m] + distances[m][j];
		/*
		 * double cost = collCost * distances[i][k] + transCost *
		 * distances[k][m] + distCost * distances[m][j];
		 */
		return cost;
	}

	/**
	 * q(k,m)
	 */
	private static double q(int k, int m) {
		if (k == m)
			return 0;
		else
			return q;
	}

	public static void main(String[] args) throws FileNotFoundException {
		// Filling in the flows matrix assymetrically
		for (int i = 0; i < nVar; i++) {
			for (int j = 0; j < nVar; j++) {
				flows[i][j] = tmpFlows[i][j] + tmpFlows[j][i];
			}
		}
		System.out.println(_R);
		try {
			GRBEnv env = new GRBEnv("RpHND_Dual.log");
			GRBModel model = new GRBModel(env);

			// Create variables
			GRBVar u2 = model.addVar(0.0, GRB.INFINITY, P, GRB.CONTINUOUS, "u2");
			GRBVar[][] u3 = new GRBVar[nVar][nVar];
			GRBVar[][][][] u4 = new GRBVar[nVar][nVar][nVar][R + 1];
			GRBVar[][][][] u5 = new GRBVar[nVar][nVar][nVar][R + 1];
			GRBVar[][][] u6 = new GRBVar[nVar][nVar][R + 1];
			GRBVar[][][] u7 = new GRBVar[nVar][nVar][R + 1];
			GRBVar[][][] u8 = new GRBVar[nVar][nVar][_R + 1]; //Leaf-nodes excluded
			GRBVar[][][] u9 = new GRBVar[nVar][nVar][_R + 1 ]; //Leaf-nodes excluded 
			GRBVar[][][][] u10 = new GRBVar[nVar][nVar][nVar][_R + 1]; //Leaf-nodes excluded 
			GRBVar[][][][] u11 = new GRBVar[nVar][nVar][nVar][_R + 1]; //Leaf-nodes excluded 

			// u3_i,j
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					String name = "u3_" + i + "_" + j;
					u3[i][j] = model.addVar(0, GRB.INFINITY, 1, GRB.CONTINUOUS, name);
				}
			}

			// u4_i,j,k,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int k = 0; k < nVar; k++) {
						for (int r = 0; r <= R; r++) {
							String name = "u4_" + i + "_" + j + "_" + k + "_" + r;
							u4[i][j][k][r] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, name);
						}
					}
				}
			}

			// u5_i,j,m,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int m = 0; m < nVar; m++) {
						for (int r = 0; r <= R; r++) {
							String name = "u5_" + i + "_" + j + "_" + m + "_" + r;
							u5[i][j][m][r] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, name);
						}
					}
				}
			}

			// u6_i,j,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= R; r++) {
						String name = "u6_" + i + "_" + j + "_" + r;
						u6[i][j][r] = model.addVar(0, GRB.INFINITY, M, GRB.CONTINUOUS, name);
					}
				}
			}

			// u7_i,j,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= R; r++) {
						String name = "u7_" + i + "_" + j + "_" + r;
						u7[i][j][r] = model.addVar(0, GRB.INFINITY, M, GRB.CONTINUOUS, name);
					}
				}
			}

			// u8_i,j,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= _R; r++) { //Leaf-nodes excluded
						String name = "u8_" + i + "_" + j + "_" + r;
						u8[i][j][r] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, name);
					}
				}
			}
			
			// u9_i,j,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= _R; r++) {	//Leaf-nodes excluded
						String name = "u9_" + i + "_" + j + "_" + r;
						u9[i][j][r] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, name);
					}
				}
			}
			
			// u10_i,j,k,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int k = 0; k < nVar; k++) {
						for (int r = 0; r <= _R; r++) { //Leaf-nodes excluded
							String name = "u10_" + i + "_" + j + "_" + k + "_" + r;
							u10[i][j][k][r] = model.addVar(0, GRB.INFINITY, M, GRB.CONTINUOUS, name);
						}
					}
				}
			}
						
			// u11_i,j,k,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int m = 0; m < nVar; m++) {
						for (int r = 0; r <= _R; r++) { //Leaf-nodes excluded
							String name = "u11_" + i + "_" + j + "_" + m + "_" + r;
							u11[i][j][m][r] = model.addVar(0, GRB.INFINITY, M, GRB.CONTINUOUS, name);
						}
					}
				}
			}
			
			// Integrate new variables
			model.update();
			
			// Set objective function
			GRBLinExpr expr = new GRBLinExpr();
			// u2
			expr.addTerm(P, u2);
			
			// u3_i,j
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					expr.addTerm(1, u3[i][j]);
				}
			}

			// u6_i,j,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= R; r++) {
						expr.addTerm(M, u6[i][j][r]);
					}
				}
			}

			// u7_i,j,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= R; r++) {
						expr.addTerm(M, u7[i][j][r]);
					}
				}
			}
			
			// u10_i,j,k,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int k = 0; k < nVar; k++) {
						for (int r = 0; r <= _R; r++) { //Leaf-nodes excluded
							expr.addTerm(M, u10[i][j][k][r]);
						}
					}
				}
			}
						
			// u11_i,j,k,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int m = 0; m < nVar; m++) {
						for (int r = 0; r <= _R; r++) { //Leaf-nodes excluded
							expr.addTerm(M, u11[i][j][m][r]);
						}
					}
				}
			}
			
			model.setObjective(expr, GRB.MAXIMIZE);
			
			// Adding constraints

			// Constraint 2			
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int k=0; k<nVar; k++){
						for (int m = 0; m < nVar; m++) {
							int r=0;
							GRBLinExpr con2 = new GRBLinExpr();
							con2.addTerm(1, u3[i][j]);
							con2.addTerm(1, u4[i][j][k][r]);
							con2.addTerm(1, u5[i][j][k][r]);
							if (i!=k && i!=m)
								con2.addTerm(1, u6[i][j][r]);
							if (j!=k && j!=m)
								con2.addTerm(1, u7[i][j][r]);
							if (r<=_R){	// if 'r' is not a leaf-node
								if (k!=i && k!=j)
									con2.addTerm(1, u8[i][j][r]);
								con2.addTerm(-1, u8[i][j][(int) Math.floor((r-1)/2)]);  // Math.floor((r-1)/2) defines parents node's index
								if (m!=i && m!=j)
									con2.addTerm(1, u9[i][j][r]);
								con2.addTerm(-1, u9[i][j][(int) Math.floor((r-1)/2)]);
								con2.addTerm(M, u10[i][j][k][r]);
								con2.addTerm(M, u11[i][j][m][r]);
							}
							double CoEf = flows[i][j] * Cikmj(i, k, m, j) * (1 - Q(i, k, m, j));
							model.addConstr(con2, GRB.LESS_EQUAL, CoEf, "c2" + i + "_" + j + "_" + k + "_" + m + "_" + r);
						}
					}
				}
			}
			
			// Constraint 3			
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int k=0; k<nVar; k++){
						for (int m = 0; m < nVar; m++) {
							for (int r=0; r<=R; r++){
								GRBLinExpr con3 = new GRBLinExpr();
								con3.addTerm(1, u4[i][j][k][r]);
								con3.addTerm(1, u5[i][j][k][r]);
								if (i!=k && i!=m)
									con3.addTerm(1, u6[i][j][r]);
								if (j!=k && j!=m)
									con3.addTerm(1, u7[i][j][r]);
								
								if (r<=_R){	// if 'r' is not a leaf-node
									if (k!=i && k!=j)
										con3.addTerm(1, u8[i][j][r]);
									
									if (m!=i && m!=j)
										con3.addTerm(1, u9[i][j][r]);
									if (r%2 == 0)
										con3.addTerm(-1, u9[i][j][(int) Math.floor((r-1)/2)]);
									else
										con3.addTerm(-1, u8[i][j][(int) Math.floor((r-1)/2)]);										
									con3.addTerm(M, u10[i][j][k][r]);
									for (int s : BinaryTree.leftChildren(r, D)){
										con3.addTerm(1, u10[i][j][k][s]);
										con3.addTerm(1, u10[i][j][m][s]);
									}
									con3.addTerm(M, u11[i][j][m][r]);
									for (int s : BinaryTree.rightChildren(r, D)){
										con3.addTerm(1, u11[i][j][k][s]);
										con3.addTerm(1, u11[i][j][m][s]);
									}
								}
								System.out.println();
								double CoEf = flows[i][j] * Cikmj(i, k, m, j) * Math.pow(q, Math.floor(Math.log(r+1)/Math.log(2)));
								model.addConstr(con3, GRB.LESS_EQUAL, CoEf, "c3" + i + "_" + j + "_" + k + "_" + m + "_" + r);
							}							
						}
					}
				}
			}
			
			// Constraint 4
			
			for (int x=0; x<nVar; x++){
				GRBLinExpr con4 = new GRBLinExpr();
				con4.addTerm(1, u2);
				for (int i = 0; i < nVar; i++) {
					for (int j = i + 1; j < nVar; j++) {
						for (int r = 0; r <= R; r++) {
							con4.addTerm(-1, u4[i][j][x][r]);
						}
					}
				}
				for (int i = 0; i < nVar; i++) {
					for (int j = i + 1; j < nVar; j++) {
						for (int r = 0; r <= R; r++) {
							con4.addTerm(-1, u5[i][j][x][r]);
						}
					}
				}
				for (int j=x+1; j<nVar; j++){
					for (int r = 0; r <= R; r++) {
						con4.addTerm(M, u6[x][j][r]);
					}
				}
				for (int i=0; i<nVar; i++){
					for (int r = 0; r <= R; r++) {
						if (x>i)
							con4.addTerm(M, u7[i][x][r]);
					}
				}
				model.addConstr(con4, GRB.LESS_EQUAL, 0, "con4" + x);
			}
			
			// Optimize model
			model.optimize();
			
			//Printing solution to a file
			File file = new File("D:/Dual_results.txt");
			PrintWriter out = new PrintWriter(file);
			
			// u3_i,j
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					out.println(u3[i][j].get(GRB.StringAttr.VarName)
	                         + " " + u3[i][j].get(GRB.DoubleAttr.X));
				}
			}

			// u4_i,j,k,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int k = 0; k < nVar; k++) {
						for (int r = 0; r <= R; r++) {
							out.println(u4[i][j][k][r].get(GRB.StringAttr.VarName)
			                         + " " + u4[i][j][k][r].get(GRB.DoubleAttr.X));
						}
					}
				}
			}

			// u5_i,j,m,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int m = 0; m < nVar; m++) {
						for (int r = 0; r <= R; r++) {
							out.println(u5[i][j][m][r].get(GRB.StringAttr.VarName)
			                         + " " + u5[i][j][m][r].get(GRB.DoubleAttr.X));
						}
					}
				}
			}

			// u6_i,j,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= R; r++) {
						out.println(u6[i][j][r].get(GRB.StringAttr.VarName)
		                         + " " + u6[i][j][r].get(GRB.DoubleAttr.X));
					}
				}
			}

			// u7_i,j,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= R; r++) {
						out.println(u7[i][j][r].get(GRB.StringAttr.VarName)
		                         + " " + u7[i][j][r].get(GRB.DoubleAttr.X));
					}
				}
			}

			// u8_i,j,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= _R; r++) { //Leaf-nodes excluded
						out.println(u8[i][j][r].get(GRB.StringAttr.VarName)
		                         + " " + u8[i][j][r].get(GRB.DoubleAttr.X));
					}
				}
			}
			
			// u9_i,j,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= _R; r++) {	//Leaf-nodes excluded
						out.println(u9[i][j][r].get(GRB.StringAttr.VarName)
		                         + " " + u9[i][j][r].get(GRB.DoubleAttr.X));
					}
				}
			}
			
			// u10_i,j,k,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int k = 0; k < nVar; k++) {
						for (int r = 0; r <= _R; r++) { //Leaf-nodes excluded
							out.println(u10[i][j][k][r].get(GRB.StringAttr.VarName)
			                         + " " + u10[i][j][k][r].get(GRB.DoubleAttr.X));
						}
					}
				}
			}
						
			// u11_i,j,k,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int m = 0; m < nVar; m++) {
						for (int r = 0; r <= _R; r++) { //Leaf-nodes excluded
							out.println(u11[i][j][m][r].get(GRB.StringAttr.VarName)
			                         + " " + u11[i][j][m][r].get(GRB.DoubleAttr.X));
						}
					}
				}
			}
			
			//Results 
			System.out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal));
			out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal));
			// Dispose of model and environment
		    model.dispose();
		    env.dispose();
			out.close();
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
		}
	}
}
