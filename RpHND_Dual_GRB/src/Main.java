import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.security.InvalidAlgorithmParameterException;
import java.util.List;

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
	private static int P = 2; // number of hubs to be located
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

	/*public static void main(String[] args) throws InvalidAlgorithmParameterException{
		
		int i = 53921;
		if (i%2 == 0)
			System.out.println("Even");
		else
			System.out.println("odd");
	}*/
	public static void main(String[] args) throws FileNotFoundException, InvalidAlgorithmParameterException {
		// Filling in the flows matrix assymetrically
		for (int i = 0; i < nVar; i++) {
			for (int j = 0; j < nVar; j++) {
				flows[i][j] = tmpFlows[i][j] + tmpFlows[j][i];
			}
		}
		
		try {
			GRBEnv env = new GRBEnv("RpHND_Dual.log");
			GRBModel model = new GRBModel(env);
			double tot = 0; 
			// Create variables
			GRBVar u2 = model.addVar(-1*GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, "u2");
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
					u3[i][j] = model.addVar(-1*GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, name);
				}
			}

			// u4_i,j,k,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int k = 0; k < nVar; k++) {
						for (int r = 0; r <= R; r++) {
							String name = "u4_" + i + "_" + j + "_" + k + "_" + r;
							u4[i][j][k][r] = model.addVar(-1*GRB.INFINITY, 0, 0, GRB.CONTINUOUS, name);
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
							u5[i][j][m][r] = model.addVar(-1*GRB.INFINITY, 0, 0, GRB.CONTINUOUS, name);
						}
					}
				}
			}

			// u6_i,j,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= R; r++) {
						String name = "u6_" + i + "_" + j + "_" + r;
						u6[i][j][r] = model.addVar(-1*GRB.INFINITY, 0, 0, GRB.CONTINUOUS, name);
					}
				}
			}

			// u7_i,j,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= R; r++) {
						String name = "u7_" + i + "_" + j + "_" + r;
						u7[i][j][r] = model.addVar(-1*GRB.INFINITY, 0, 0, GRB.CONTINUOUS, name);
					}
				}
			}

			// u8_i,j,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= _R; r++) { //Leaf-nodes excluded
						String name = "u8_" + i + "_" + j + "_" + r;
						u8[i][j][r] = model.addVar(-1*GRB.INFINITY, 0, 0, GRB.CONTINUOUS, name);
					}
				}
			}
			
			// u9_i,j,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= _R; r++) {	//Leaf-nodes excluded
						String name = "u9_" + i + "_" + j + "_" + r;
						u9[i][j][r] = model.addVar(-1*GRB.INFINITY, 0, 0, GRB.CONTINUOUS, name);
					}
				}
			}
			
			// u10_i,j,k,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int k = 0; k < nVar; k++) {
						for (int r = 0; r <= _R; r++) { //Leaf-nodes excluded
							String name = "u10_" + i + "_" + j + "_" + k + "_" + r;
							u10[i][j][k][r] = model.addVar(-1*GRB.INFINITY, 0, 0, GRB.CONTINUOUS, name);
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
							u11[i][j][m][r] = model.addVar(-1*GRB.INFINITY, 0, 0, GRB.CONTINUOUS, name);
						}
					}
				}
			}
			
			// Integrate new variables
			model.update();
			
			// Set objective function
			GRBLinExpr expr = new GRBLinExpr();
			PrintWriter out = new PrintWriter(new File("D:/modelDual.txt"));
			out.println("Obj");
			// u2
			expr.addTerm(P, u2);
			out.println("+ " + P + " " + u2.get(GRB.StringAttr.VarName));
			
			// u3_i,j
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					expr.addTerm(1, u3[i][j]);
					out.println("+ " + 1 + " " + u3[i][j].get(GRB.StringAttr.VarName));
				}
			}

			// u6_i,j,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= R; r++) {
						expr.addTerm(M, u6[i][j][r]);
						out.println("+ " + M + " " + u6[i][j][r].get(GRB.StringAttr.VarName));
					}
				}
			}

			// u7_i,j,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= R; r++) {
						expr.addTerm(M, u7[i][j][r]);
						out.println("+ " + M + " " + u7[i][j][r].get(GRB.StringAttr.VarName));
					}
				}
			}
			
			// u10_i,j,k,r
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int k = 0; k < nVar; k++) {
						for (int r = 0; r <= _R; r++) { //Leaf-nodes excluded
							expr.addTerm(M, u10[i][j][k][r]);
							out.println("+ " + M + " " + u10[i][j][k][r].get(GRB.StringAttr.VarName));
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
							out.println("+ " + M + " " + u11[i][j][m][r].get(GRB.StringAttr.VarName));
						}
					}
				}
			}
			
			model.setObjective(expr, GRB.MAXIMIZE);
			
			// Adding constraints

			// Constraint 2	
			out.println("Constraint2");
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int k=0; k<nVar; k++){
						for (int m = 0; m < nVar; m++) {
							int r=0;
							GRBLinExpr con2 = new GRBLinExpr();
							con2.addTerm(1, u3[i][j]); out.print(" + " + 1 + u3[i][j].get(GRB.StringAttr.VarName));
							con2.addTerm(1, u4[i][j][k][r]); out.print(" + " + 1 + u4[i][j][k][r].get(GRB.StringAttr.VarName));
							con2.addTerm(1, u5[i][j][k][r]); out.print(" + " + 1 + u5[i][j][k][r].get(GRB.StringAttr.VarName));
							if (i!=k && i!=m)
								{tot+=1;con2.addTerm(1, u6[i][j][r]); out.print(" + " + 1 + u6[i][j][r].get(GRB.StringAttr.VarName));}
							if (j!=k && j!=m)
								{tot+=1;con2.addTerm(1, u7[i][j][r]); out.print(" + " + 1 + u7[i][j][r].get(GRB.StringAttr.VarName));}
							if (k!=i && k!=j)
								{tot+=1;con2.addTerm(1, u8[i][j][r]); out.print(" + " + 1 + u8[i][j][r].get(GRB.StringAttr.VarName));}
//							con2.addTerm(-1, u8[i][j][(int) Math.floor((r-1)/2)]); out.print(" - " + 1 + u8[i][j][(int) Math.floor((r-1)/2)].get(GRB.StringAttr.VarName)); // Math.floor((r-1)/2) defines parents node's index
							if (m!=i && m!=j)
								{tot+=1;con2.addTerm(1, u9[i][j][r]); out.print(" + " + 1 + u9[i][j][r].get(GRB.StringAttr.VarName));}
//							con2.addTerm(-1, u9[i][j][(int) Math.floor((r-1)/2)]); out.print(" - " + 1 + u9[i][j][(int) Math.floor((r-1)/2)].get(GRB.StringAttr.VarName));
							con2.addTerm(M, u10[i][j][k][r]); out.print(" + " + M + u10[i][j][k][r].get(GRB.StringAttr.VarName));
							con2.addTerm(M, u11[i][j][m][r]); out.print(" + " + M + u11[i][j][m][r].get(GRB.StringAttr.VarName));
							double CoEf = flows[i][j] * Cikmj(i, k, m, j) * (1 - Q(i, k, m, j));
							model.addConstr(con2, GRB.LESS_EQUAL, CoEf, "c2" + i + "_" + j + "_" + k + "_" + m + "_" + r); out.println(" <= " + CoEf + "  " + i + "_" + k + "_" + m + "_" + j + "_" + r);
							tot +=(1+1+1+M+M);
						}
					}
				}
			}
			
			// Constraint 3		
			out.println("Constraint3");
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int k=0; k<nVar; k++){
						for (int m = 0; m < nVar; m++) {
							for (int r=1; r<=_R; r++){  // only non-leaf nodes are included in these constraints
								GRBLinExpr con3 = new GRBLinExpr();
								con3.addTerm(1, u4[i][j][k][r]); out.print(" + " + 1 + u4[i][j][k][r].get(GRB.StringAttr.VarName));
								con3.addTerm(1, u5[i][j][m][r]); out.print(" + " + 1 + u5[i][j][m][r].get(GRB.StringAttr.VarName));
								if (i!=k && i!=m)
									{tot+=1;con3.addTerm(1, u6[i][j][r]); out.print(" + " + 1 + u6[i][j][r].get(GRB.StringAttr.VarName));}
								if (j!=k && j!=m)
									{tot+=1;con3.addTerm(1, u7[i][j][r]); out.print(" + " + 1 +u7[i][j][r].get(GRB.StringAttr.VarName));}								
									
								if (k!=i && k!=j)
									{tot+=1;con3.addTerm(1, u8[i][j][r]); out.print(" + " + 1 + u8[i][j][r].get(GRB.StringAttr.VarName));}
								if (m!=i && m!=j)
									{tot+=1;con3.addTerm(1, u9[i][j][r]); out.print(" + " + 1 + u9[i][j][r].get(GRB.StringAttr.VarName));}
								if (r%2 == 0)
									{tot+=-1;con3.addTerm(-1, u9[i][j][(int) Math.floor((r-1)/2)]); out.print(" - " + 1 + u9[i][j][(int) Math.floor((r-1)/2)].get(GRB.StringAttr.VarName));}
								else
									{tot+=-1;con3.addTerm(-1, u8[i][j][(int) Math.floor((r-1)/2)]); out.print(" - " + 1 + u8[i][j][(int) Math.floor((r-1)/2)].get(GRB.StringAttr.VarName));}										
								con3.addTerm(M, u10[i][j][k][r]); out.print(" + " + M + u10[i][j][k][r].get(GRB.StringAttr.VarName));
								con3.addTerm(M, u11[i][j][m][r]); out.print(" + " + M + u11[i][j][m][r].get(GRB.StringAttr.VarName));
								tot+=(2+M+M);
								for (int parent : BinaryTree.getParents(r)){
									if (BinaryTree.isLeftChild(parent, r))
										{tot+=1;con3.addTerm(1, u10[i][j][k][parent]); out.print(" + " + 1 + u10[i][j][k][parent].get(GRB.StringAttr.VarName));
										tot+=1;con3.addTerm(1, u10[i][j][m][parent]); out.print(" + " + 1 + u10[i][j][m][parent].get(GRB.StringAttr.VarName));}
									else
										{tot+=1;con3.addTerm(1, u11[i][j][k][parent]); out.print(" + " + 1 + u11[i][j][k][parent].get(GRB.StringAttr.VarName));
										tot+=1;con3.addTerm(1, u11[i][j][m][parent]); out.print(" + " + 1 + u11[i][j][m][parent].get(GRB.StringAttr.VarName));}
								}
																
								System.out.println();
								double CoEf = flows[i][j] * Cikmj(i, k, m, j) * Math.pow(q, Math.floor(Math.log(r+1)/Math.log(2)));
								model.addConstr(con3, GRB.LESS_EQUAL, CoEf, "c3" + i + "_" + j + "_" + k + "_" + m + "_" + r); out.println(" <= " + CoEf + "  " + i + "_" + k + "_" + m + "_" + j + "_" + r);
							}							
						}
					}
				}
			}
			
			// Constraint 4
			out.println("constraint4");
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int k=0; k<nVar; k++){
						for (int m = 0; m < nVar; m++) {
							for (int r=_R+1; r<=R; r++){  // only leaf nodes are included in these constraints
								GRBLinExpr con4 = new GRBLinExpr();
								tot+=1;con4.addTerm(1, u4[i][j][k][r]); out.print(" + " + 1 + u4[i][j][k][r].get(GRB.StringAttr.VarName));
								tot+=1;con4.addTerm(1, u5[i][j][k][r]); out.print(" + " + 1 + u5[i][j][k][r].get(GRB.StringAttr.VarName));
								if (i!=k && i!=m)
									{tot+=1;con4.addTerm(1, u6[i][j][r]); out.print(" + " + 1 + u6[i][j][r].get(GRB.StringAttr.VarName));}
								if (j!=k && j!=m)
									{tot+=1;con4.addTerm(1, u7[i][j][r]); out.print(" + " + 1 + u7[i][j][r].get(GRB.StringAttr.VarName));}
							if (r%2 == 0)
								{tot+=-1;con4.addTerm(-1, u9[i][j][(int) Math.floor((r-1)/2)]); out.print(" - " + 1 + u9[i][j][(int) Math.floor((r-1)/2)].get(GRB.StringAttr.VarName));}
							else
								{tot+=-1;con4.addTerm(-1, u8[i][j][(int) Math.floor((r-1)/2)]); out.print(" - " + 1 + u8[i][j][(int) Math.floor((r-1)/2)].get(GRB.StringAttr.VarName));}										
		
							for (int parent : BinaryTree.getParents(r)){
								if (BinaryTree.isLeftChild(parent, r))
									{tot+=1;con4.addTerm(1, u10[i][j][k][parent]); out.print(" + " + 1 + u10[i][j][k][parent].get(GRB.StringAttr.VarName));
									tot+=1;con4.addTerm(1, u10[i][j][m][parent]); out.print(" + " + 1 + u10[i][j][m][parent].get(GRB.StringAttr.VarName));}
								else
									{tot+=1;con4.addTerm(1, u11[i][j][k][parent]); out.print(" + " + 1 + u11[i][j][k][parent].get(GRB.StringAttr.VarName));
									tot+=1;con4.addTerm(1, u11[i][j][m][parent]); out.print(" + " + 1 + u11[i][j][m][parent].get(GRB.StringAttr.VarName));}
							}
								double CoEf = flows[i][j] * Cikmj(i, k, m, j) * Math.pow(q, Math.floor(Math.log(r+1)/Math.log(2)));
								model.addConstr(con4, GRB.LESS_EQUAL, CoEf, "c3" + i + "_" + j + "_" + k + "_" + m + "_" + r); out.println(" <= " + CoEf + "  " + i + "_" + k + "_" + m + "_" + j + "_" + r);
							}							
						}
					}
				}
			}
			
			// Constraint 5
			out.println("constraint5");
			for (int x=0; x<nVar; x++){
				GRBLinExpr con5 = new GRBLinExpr();
				con5.addTerm(1, u2); out.print(" + " + 1 + u2.get(GRB.StringAttr.VarName));
				for (int i = 0; i < nVar; i++) {
					for (int j = i + 1; j < nVar; j++) {
						for (int r = 0; r <= R; r++) {
							tot+=-1;con5.addTerm(-1, u4[i][j][x][r]); out.print(" - " + 1 + u4[i][j][x][r].get(GRB.StringAttr.VarName));
						}
					}
				}
				for (int i = 0; i < nVar; i++) {
					for (int j = i + 1; j < nVar; j++) {
						for (int r = 0; r <= R; r++) {
							tot+=-1;con5.addTerm(-1, u5[i][j][x][r]); out.print(" - " + 1 + u5[i][j][x][r].get(GRB.StringAttr.VarName));
						}
					}
				}
				for (int j=x+1; j<nVar; j++){
					for (int r = 0; r <= R; r++) {
						tot+=M;con5.addTerm(M, u6[x][j][r]);  out.print(" + " + M + u6[x][j][r].get(GRB.StringAttr.VarName));
					}
				}
				for (int i=0; i<nVar; i++){
					for (int r = 0; r <= R; r++) {
						if (x>i)
							{tot+=M;con5.addTerm(M, u7[i][x][r]); out.print(" + " + M + u7[i][x][r].get(GRB.StringAttr.VarName));}
					}
				}
				model.addConstr(con5, GRB.LESS_EQUAL, 0, "con4" + x); out.println(" <= 0");
			}
			
			// Optimize model
			model.optimize();
			
			System.out.println("Number of varibles: " + model.get(GRB.IntAttr.NumVars));
			System.out.println("Number of constraints: " + model.get(GRB.IntAttr.NumConstrs));
						
			//Printing solution to a file
			File file = new File("D:/Dual_results.txt");
			PrintWriter out2 = new PrintWriter(file);
			// Dispose of model and environment
			GRBVar[] vars = model.getVars();
			for (GRBVar var: vars){
				out2.println(var.get(GRB.StringAttr.VarName) + " " + var.get(GRB.DoubleAttr.X));
			}
			out2.close();
		    model.dispose();
		    env.dispose();
			out.close();
			System.out.println(tot);
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
		}
	}
}
