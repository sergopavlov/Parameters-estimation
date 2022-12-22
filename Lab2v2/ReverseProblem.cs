using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Lab2v2
{
    public class ReverseProblem
    {
        public List<Parameter> parameters;
        public List<double> TrueComponentValues;
        public DirectProblem problem;
        public int ComponentCount;
        public double alpha { get; set; }
        public double CalcF()
        {
            double res = 0;
            for (int i = 0; i < ComponentCount; i++)
            {
                res += (problem.CalcValueAtPoint(i, parameters) - TrueComponentValues[i]) * (problem.CalcValueAtPoint(i, parameters) - TrueComponentValues[i]) / (TrueComponentValues[i] * TrueComponentValues[i]);
            };
            return res;
        }
        public void Solve(int maxiter, double eps)
        {
            List<double> weights = new();
            for (int i = 0; i < ComponentCount; i++)
            {
                weights.Add(1 / TrueComponentValues[i]);
            }
            MatrixFull mat = new MatrixFull(parameters.Count);
            List<double> b = new();
            for (int i = 0; i < parameters.Count; i++)
            {
                b.Add(0);
            }
            double penalty = CalcF();
            int k = 0;
            while (penalty > eps && k < maxiter)
            {
                List<List<double>> diffs = problem.DiffComponentByParameter(parameters);
                double[] V0 = new double[ComponentCount];
                for (int i = 0; i < ComponentCount; i++)
                {
                    V0[i] = problem.CalcValueAtPoint(i, parameters);
                }

                for (int s = 0; s < ComponentCount; s++)
                {
                    for (int i = 0; i < parameters.Count; i++)
                    {
                        b[i] += weights[s] * weights[s] * diffs[s][i] * (TrueComponentValues[s] - V0[s]);
                        for (int j = 0; j < parameters.Count; j++)
                        {
                            mat.mat[i][j] += weights[s] * weights[s] * diffs[s][i] * diffs[s][j];
                        }
                    }
                }
                for (int i = 0; i < parameters.Count; i++)
                {
                    mat.mat[i][i] += alpha;
                }
                //добавить регуляризацию
                var dparam = mat.SolveLU(b);
                mat.Clear();
                for (int i = 0; i < parameters.Count; i++)
                {
                    parameters[i].Value += dparam[i];
                    b[i] = 0;
                }
                penalty = CalcF();
                Console.WriteLine($"{penalty}";
                foreach (var item in parameters)
                {
                    Console.WriteLine(item);
                }
                k++;
            }
        }
    }
}
