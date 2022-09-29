using System;

namespace myns
{
    internal class Program
    {
        static void Main(string[] args)
        {
            List<Vector> A = new();
            List<Vector> B = new();
            List<Vector> M = new();
            List<Vector> N = new();
            List<double> V = new();
            foreach (var item in File.ReadLines("txt/напряжения.txt"))
            {
                V.Add(double.Parse(item));
            }
            foreach (var item in File.ReadAllLines("txt/источники.txt"))
            {
                var cur = item.Split(' ');
                A.Add(new Vector(double.Parse(cur[0]), double.Parse(cur[1])));
                B.Add(new Vector(double.Parse(cur[2]), double.Parse(cur[3])));
            }
            foreach (var item in File.ReadAllLines("txt/приемники.txt"))
            {
                var cur = item.Split(' ');
                M.Add(new Vector(double.Parse(cur[0]), double.Parse(cur[1])));
                N.Add(new Vector(double.Parse(cur[2]), double.Parse(cur[3])));
            }
            double sigma = double.Parse(File.ReadAllText("txt/sigma.txt"));
            int n = A.Count;
            int m = M.Count;
            Matrix mat = new Matrix(n);
            List<double> b = new();
            List<double> res = new();
            for (int i = 0; i < n; i++)
            {
                b.Add(0);
                res.Add(0);
            }
            
            double penalty = 0;
            for (int i = 0; i < m; i++)
            {
                penalty += (CalcV(i,res,sigma,A,B,M,N) - V[i])* (CalcV(i, res, sigma, A, B, M, N) - V[i]) / (V[i] * V[i]);
            }
            penalty = Math.Sqrt(penalty);
            Console.WriteLine($"{res[0]} {res[1]} {res[2]} {penalty}");
            
            while (penalty > 1e-14)
            {


                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        for (int k = 0; k < m; k++)
                        {
                            mat.mat[i][j] += 1 / (V[k] * V[k] * 4 * Math.PI * Math.PI * sigma * sigma) * (1 / (B[j] - M[k]).Norm - 1 / (A[j] - M[k]).Norm - 1 / (B[j] - N[k]).Norm + 1 / (A[j] - N[k]).Norm) * (1 / (B[i] - M[k]).Norm - 1 / (A[i] - M[k]).Norm - 1 / (B[i] - N[k]).Norm + 1 / (A[i] - N[k]).Norm);
                        }
                    }
                    b[i] = 0;
                    for (int k = 0; k < m; k++)
                    {
                        b[i] += 1 / (V[k] * V[k] * 2 * Math.PI * sigma) * (1 / (B[i] - M[k]).Norm - 1 / (A[i] - M[k]).Norm - 1 / (B[i] - N[k]).Norm + 1 / (A[i] - N[k]).Norm) * (V[k] - CalcV(k, res, sigma, A, B, M, N));
                    }
                    mat.mat[i][i] += 1e-15;
                }

                var dres= mat.SolveLU(b);
                for (int i = 0; i < n; i++)
                {
                    res[i] += dres[i];
                }
                mat.Clear();
                penalty = 0;
                for (int i = 0; i < m; i++)
                {
                    penalty += (CalcV(i, res, sigma, A, B, M, N) - V[i]) * (CalcV(i, res, sigma, A, B, M, N) - V[i]) / (V[i] * V[i]);
                }
                penalty = Math.Sqrt(penalty);
                Console.WriteLine($"{res[0]} {res[1]} {res[2]} {penalty}");
            }
           

        }

        static double CalcV(int i, List<double> I, double sigma, List<Vector> A, List<Vector> B, List<Vector> M, List<Vector> N)
        {
            double res = 0;
            int n = I.Count;
            for (int j = 0; j < n; j++)
            {
                res += 1 / (2 * Math.PI * sigma) * I[j] * (1 / (B[j] - M[i]).Norm - 1 / (A[j] - M[i]).Norm - 1 / (B[j] - N[i]).Norm + 1 / (A[j] - N[i]).Norm);
            }
            return res;
        }
    }

    public class Matrix
    {

        public void Clear()
        {
            for (int i = 0; i < dim; i++)
            {
                for (int j = 0; j < dim; j++)
                {
                    mat[i][j] = 0;
                }
            }
        }
        public List<double> Mult(List<double> x)
        {
            List<double> res = new();
            for (int i = 0; i < dim; i++)
            {
                res.Add(0);
            }
            for (int i = 0; i < dim; i++)
            {
                for (int j = 0; j < dim; j++)
                {
                    res[i] += mat[i][j] * x[j];
                }
            }
            return res;
        }
        public List<List<double>> mat = new();
        public int dim { get; init; }
        public Matrix(int n)
        {
            dim = n;
            for (int i = 0; i < n; i++)
            {
                mat.Add(new List<double>());
                for (int j = 0; j < n; j++)
                {
                    mat[i].Add(0);
                }
            }
        }
        void LUDecomposition()
        {
            int n = dim;
            for (int i = 0; i < n; i++)
            {
                double summd = 0;
                for (int j = 0; j < i; j++)
                {
                    double summl = 0;
                    double summu = 0;
                    for (int k = 0; k < j; k++)
                    {
                        summl += mat[i][k] * mat[k][j];
                        summu += mat[j][k] * mat[k][i];
                    }
                    mat[i][j] -= summl;
                    mat[j][i] = (mat[j][i] - summu) / mat[j][j];
                    summd += mat[i][j] * mat[j][i];
                }
                mat[i][i] -= summd;
            }
        }
        void Gauss(List<double> b)
        {
            int n = dim;
            for (int i = 0; i < n; i++)
            {
                double summ = 0;
                for (int j = 0; j < i; j++)
                {
                    summ += mat[i][j] * b[j];
                }
                b[i] = (b[i] - summ) / mat[i][i];
            }
            for (int i = n - 1; i >= 0; i--)
            {
                double summ = 0;
                for (int j = n - 1; j > i; j--)
                {
                    summ += mat[i][j] * b[j];
                }
                b[i] -= summ;
            }
        }
        public List<double> SolveLU(List<double> b)
        {
            LUDecomposition();
            Gauss(b);
            return b;
        }
    }

    public class Vector
    {
        public Vector(double x, double y)
        {
            this.x = x;
            this.y = y;
        }

        public static Vector operator +(Vector v1, Vector v2) => new Vector(v1.x + v2.x, v1.y + v2.y);
        public static Vector operator -(Vector v1, Vector v2) => new Vector(v1.x - v2.x, v1.y - v2.y);
        public static double operator *(Vector v1, Vector v2) => v1.x * v2.x + v1.y * v2.y;
        public double x { get; init; }
        public double y { get; init; }
        public double Norm
        {
            get => Math.Sqrt(x * x + y * y);
        }
    }
}