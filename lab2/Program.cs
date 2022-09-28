using System;
using System.IO;

namespace lab2
{
    internal class Program
    {
        static void Main(string[] args)
        {
            MFE mfe = new MFE("txt/Grid.txt", 30, 1, 1);
            Console.WriteLine($"{mfe.GetSollution(1, 0)} {mfe.GetSollution(2, 0)} {mfe.GetSollution(5, 0)} {mfe.GetSollution(10, 0)}");
            Console.WriteLine($"{mfe.GetSollution(1, 0)} {mfe.GetSollution(2, 0)} {mfe.GetSollution(5, 0)} {mfe.GetSollution(10, 0)}");
        }

    }
    public class MFE
    {
        public MFE(string Grid, double h, double sigma1, double sigma2)
        {
            this.h = h;
            this.sigma1 = sigma1;
            this.sigma2 = sigma2;
            var text = File.ReadAllLines(Grid);
            var cur = text[0].Split(' ');
            var r0 = double.Parse(cur[0]);
            var r1 = double.Parse(cur[1]);
            rn = int.Parse(cur[2]);
            var rq = double.Parse(cur[3]);
            var rl = (r1 - r0) * (1 - rq) / (1 - Math.Pow(rq, rn));
            if (rq == 1)
                rl = (r1 - r0) / rn;
            rgrid.Add(r0);
            for (int i = 0; i < rn; i++, rl *= rq)
            {
                rgrid.Add(rgrid[i] + rl);
            }
            cur = text[1].Split(' ');
            var z0 = double.Parse(cur[0]);
            var z1 = double.Parse(cur[1]);
            zn = int.Parse(cur[2]);
            var zq = double.Parse(cur[3]);
            var zl = (z1 - z0) * (1 - zq) / (1 - Math.Pow(zq, zn));
            if (zq == 1)
                zl = (z1 - z0) / zn;
            zgrid.Add(z0);
            for (int i = 0; i < zn; i++, zl *= zq)
            {
                zgrid.Add(zgrid[i] + zl);
            }
            for (int i = 0; i <= zn; i++)
            {
                zgrid[i] = -rgrid[zn - i];
            }
            mat = new Matrix((rn + 1) * (zn + 1));
            Console.WriteLine($"{rgrid[1] - rgrid[0]} {rgrid[100] - rgrid[99]}");
            GenerateProfile(rn, zn);
        }
        private void AddElems()
        {
            for (int r = 0; r < rn; r++)
            {
                for (int z = 0; z < zn; z++)
                {
                    Dictionary<int, int> dic = new()
                    {
                        { 0, r + (rn + 1) * z },
                        { 1, r + 1 + (rn + 1) * z },
                        { 2, r + (rn + 1) * (z + 1) },
                        { 3, r + 1 + (rn + 1) * (z + 1) },
                    };
                    double hr = rgrid[r + 1] - rgrid[r];
                    double hz = zgrid[z + 1] - zgrid[z];
                    double cursigma = (zgrid[z + 1] + zgrid[z]) / 2 > h ? sigma1 : sigma2;
                    double[] rs = new double[] { rgrid[r], rgrid[r + 1] };
                    for (int k = 0; k < 2; k++)
                    {
                        for (int i = 0; i < 4; i++)
                        {
                            mat.di[dic[i]] += cursigma * rs[k] * (hz / hr * Matrices.GR[k][i % 2][i % 2] * Matrices.MZ[i / 2][i / 2] + hr / hz * Matrices.MR[k][i % 2][i % 2] * Matrices.GZ[i / 2][i / 2]);
                            for (int j = 0; j < i; j++)
                            {
                                //mat.b[dic[i]] += hz * hr * Matrices.MR[i%2][j%2][k]*Matrices.MZ[i/2][j/2];
                                int index = mat.ia[dic[i]];
                                while (mat.ja[index] < dic[j])
                                    index++;
                                mat.al[index] += cursigma * rs[k] * (hz / hr * Matrices.GR[k][i % 2][j % 2] * Matrices.MZ[i / 2][j / 2] + hr / hz * Matrices.MR[k][i % 2][j % 2] * Matrices.GZ[i / 2][j / 2]);
                            }
                        }
                    }
                    if (z == zn - 1 && r == 0)
                    {
                        mat.b[r + (rn + 1) * (z + 1)] += 1;
                    }
                }
            }
        }
        private void AddBoundary()
        {


            for (int r = 0; r <= rn; r++)
            {
                mat.di[r] = 1;
                mat.b[r] = 0;
                for (int k = mat.ia[r]; k < mat.ia[r + 1]; k++)
                {
                    //mat.b[mat.ja[k]] -= mat.al[k];
                    mat.al[k] = 0;
                }
                for (int i = r + 1; i < mat.n; i++)
                {
                    int k = mat.ia[i];
                    while (mat.ja[k] < r && k < mat.ia[i + 1] - 1)
                    {
                        k++;
                    }
                    if (mat.ja[k] == r)
                    {
                        //mat.b[i] -= mat.al[k];
                        mat.al[k] = 0;
                    }
                }
            }
            for (int z = 0; z <= zn; z++)
            {
                int index = rn + (rn + 1) * z;
                mat.di[index] = 1;
                mat.b[index] = 0;
                for (int k = mat.ia[index]; k < mat.ia[index + 1]; k++)
                {
                    //mat.b[mat.ja[k]] -= mat.al[k];
                    mat.al[k] = 0;
                }
                for (int i = index + 1; i < mat.n; i++)
                {
                    int k = mat.ia[i];
                    while (mat.ja[k] < index && k < mat.ia[i + 1] - 1)
                    {
                        k++;
                    }
                    if (mat.ja[k] == index)
                    {
                        //mat.b[i] -= mat.al[k];
                        mat.al[k] = 0;
                    }
                }
            }
        }
        private int rn;
        private int zn;
        private static class Matrices
        {
            public static double[][][] GR = new double[][][]
            {
                new double[][]
                {
                    new double[]{0.5,-0.5},
                    new double[]{-0.5,0.5}
                },
                new double[][]
                {
                    new double[]{0.5,-0.5},
                    new double[]{-0.5,0.5}
                }
            };
            public static double[][][] MR = new double[][][]
            {
                new double[][]
                {
                    new double[]{0.25,1.0/12},
                    new double[]{ 1.0 / 12, 1.0 / 12 }
                },
                new double[][]
                {
                    new double[]{ 1.0 / 12, 1.0/12},
                    new double[]{ 1.0 / 12, 0.25}
                }
            };
            public static double[][] MZ = new double[][]
            {
                new double[]{1.0/3,1.0/6},
                new double[]{ 1.0/6, 1.0 / 3 }
            };
            public static double[][] GZ = new double[][]
            {
                new double[]{1,-1},
                new double[]{ -1,1}
            };
        }
        public bool IsSolved { get; private set; }
        private void GenerateProfile(int n, int m)
        {
            int dim = mat.n;
            List<SortedSet<int>> list = new List<SortedSet<int>>(dim);
            for (int i = 0; i < dim; i++)
            {
                list.Add(new SortedSet<int>());
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    list[i + (n + 1) * j].Add(i + (n + 1) * j);
                    list[i + (n + 1) * j].Add(i + (n + 1) * (j + 1));
                    list[i + (n + 1) * j].Add((i + 1) + (n + 1) * j);
                    list[i + (n + 1) * j].Add((i + 1) + (n + 1) * (j + 1));

                    list[i + (n + 1) * (j + 1)].Add(i + (n + 1) * j);
                    list[i + (n + 1) * (j + 1)].Add(i + (n + 1) * (j + 1));
                    list[i + (n + 1) * (j + 1)].Add((i + 1) + (n + 1) * j);
                    list[i + (n + 1) * (j + 1)].Add((i + 1) + (n + 1) * (j + 1));

                    list[i + 1 + (n + 1) * j].Add(i + (n + 1) * j);
                    list[i + 1 + (n + 1) * j].Add(i + (n + 1) * (j + 1));
                    list[i + 1 + (n + 1) * j].Add((i + 1) + (n + 1) * j);
                    list[i + 1 + (n + 1) * j].Add((i + 1) + (n + 1) * (j + 1));

                    list[i + 1 + (n + 1) * (j + 1)].Add(i + (n + 1) * j);
                    list[i + 1 + (n + 1) * (j + 1)].Add(i + (n + 1) * (j + 1));
                    list[i + 1 + (n + 1) * (j + 1)].Add((i + 1) + (n + 1) * j);
                    list[i + 1 + (n + 1) * (j + 1)].Add((i + 1) + (n + 1) * (j + 1));
                }
            }
            mat.ia.Add(0);
            mat.ia.Add(0);
            mat.di.Add(0);
            mat.b.Add(0);
            for (int i = 1; i < dim; i++)
            {
                mat.di.Add(0);
                mat.b.Add(0);
                int count = 0;
                foreach (var item in list[i])
                {
                    if (item < i)
                    {
                        mat.ja.Add(item);
                        mat.al.Add(0);
                        count++;
                    }
                }
                mat.ia.Add(mat.ia[i] + count);
            }
        }
        public double GetSollution(double r, double z)
        {
            if (!IsSolved)
                Solve();
            if (r < 0)
                throw new Exception();
            int indexr = 1;
            while (rgrid[indexr] < r)
                indexr++;
            indexr--;
            int indexz = 1;
            while (zgrid[indexz] < z)
                indexz++;
            indexz--;
            Dictionary<int, int> dic = new()
            {
                { 0, indexr + (rn + 1) * indexz },
                { 1, indexr + 1 + (rn + 1) * indexz },
                { 2, indexr + (rn + 1) * (indexz + 1) },
                { 3, indexr + 1 + (rn + 1) * (indexz + 1) },
            };
            double rloc = (r - rgrid[indexr]) / (rgrid[indexr + 1] - rgrid[indexr]);
            double zloc = (z - zgrid[indexz]) / (zgrid[indexz + 1] - zgrid[indexz]);
            return q[dic[0]] * rloc * zloc + q[dic[1]] * (1 - rloc) * zloc + q[dic[2]] * rloc * (1 - zloc) + q[dic[3]] * (1 - rloc) * (1 - zloc);
        }
        private double h;
        public double H
        {
            get => h;
            set
            {
                h = value;
                IsSolved = false;
            }
        }
        private Matrix mat;
        private double sigma1;
        public double Sigma1
        {
            get => sigma1;
            set
            {
                sigma1 = value;
                IsSolved = false;
            }
        }
        private double sigma2;
        public double Sigma2
        {
            get => sigma2;
            set
            {
                sigma2 = value;
                IsSolved = false;
            }
        }
        private List<double> q;
        private List<double> rgrid = new();
        private List<double> zgrid = new();
        public void Solve()
        {
            if (IsSolved)
                return;
            mat.Clear();
            AddElems();
            AddBoundary();
            q = mat.SolveSLAE();
            IsSolved = true;
            return;
        }
    }
    public class Matrix
    {

        public int n;
        public List<int> ia = new();
        public List<int> ja = new();
        public List<double> di = new();
        public List<double> al = new();
        public List<double> b = new();
        public Matrix(int n)
        {
            this.n = n;
        }
        public void Clear()
        {
            for (int i = 0; i < n; i++)
            {
                di[i] = 0;
                b[i] = 0;
            }
            int m = al.Count;
            for (int i = 0; i < m; i++)
            {
                al[i] = 0;
            }
        }
        private List<double> MSG(List<double> x0, int maxiter, double eps)
        {
            double bnorm = Math.Sqrt(DotProduct(b, b));
            List<double> p = new List<double>(n);
            List<double> q = new List<double>(n);
            for (int i = 0; i < n; i++)
            {
                p.Add(0);
            }
            var r = MatrixMult(x0);
            for (int i = 0; i < n; i++)
            {
                r[i] = b[i] - r[i];
                p[i] = r[i];
            }
            int k = 0;
            double alpha, betta, rnorm = Math.Sqrt(DotProduct(r, r));
            while (k < maxiter && rnorm / bnorm > eps)
            {
                q = MatrixMult(p);
                alpha = DotProduct(r, r) / DotProduct(q, p);
                betta = 1 / DotProduct(r, r);
                for (int i = 0; i < n; i++)
                {
                    x0[i] += alpha * p[i];
                    r[i] -= alpha * q[i];
                }
                rnorm = Math.Sqrt(DotProduct(r, r));
                betta *= DotProduct(r, r);
                for (int i = 0; i < n; i++)
                {
                    p[i] = r[i] + betta * p[i];
                }
                k++;
            }
            Console.WriteLine($"{k} {rnorm / bnorm}");
            return x0;
        }
        public List<double> LOS(List<double> x0, int maxiter, double eps)
        {
            double bnorm = Math.Sqrt(DotProduct(b, b));

            List<double> Ar = new();
            List<double> r = MatrixMult(x0);
            List<double> z = new();
            for (int i = 0; i < n; i++)
            {
                r[i] = b[i] - r[i];
                z.Add(r[i]);
            }
            List<double> p = MatrixMult(z);
            int k = 0;
            double alpha, betta, rnorm = Math.Sqrt(DotProduct(r, r));
            while (k < maxiter && rnorm / bnorm > eps)
            {
                alpha = DotProduct(p, r) / DotProduct(p, p);
                for (int i = 0; i < n; i++)
                {
                    x0[i] += alpha * z[i];
                    r[i] -= alpha * p[i];
                }
                Ar = MatrixMult(r);
                betta = -DotProduct(p, Ar) / DotProduct(p, p);
                rnorm = Math.Sqrt(DotProduct(r, r));
                for (int i = 0; i < n; i++)
                {
                    z[i] = r[i] + betta * z[i];
                    p[i] = Ar[i] + betta * p[i];
                }
                k++;
            }
            Console.WriteLine($"{rnorm / bnorm} {k}");
            return x0;
        }
        double DotProduct(List<double> vec1, List<double> vec2)
        {
            if (vec1.Count != vec2.Count)
                throw new Exception();
            double res = 0;
            int m = vec1.Count;
            for (int i = 0; i < m; i++)
            {
                res += vec1[i] * vec2[i];
            }
            return res;
        }
        List<double> MatrixMult(List<double> x)
        {
            if (x.Count != n)
                throw new Exception();
            List<double> res = new List<double>(x.Count);
            for (int i = 0; i < n; i++)
            {
                res.Add(0);
            }
            for (int i = 0; i < n; i++)
            {
                res[i] = x[i] * di[i];
                for (int k = ia[i]; k < ia[i + 1]; k++)
                {
                    int j = ja[k];
                    res[i] += al[k] * x[j];
                    res[j] += al[k] * x[i];
                }
            }
            return res;
        }
        public List<double> SolveSLAE()
        {
            List<double> x0 = new();
            for (int i = 0; i < n; i++)
            {
                x0.Add(0);
            }
            return LOS(x0, 10000, 1e-15);
        }
        public List<double> SolveSLAE(List<double> x0)
        {
            for (int i = 0; i < n; i++)
            {
                x0.Add(0);
            }
            return MSG(x0, 10000, 1e-15);
        }
    }
}