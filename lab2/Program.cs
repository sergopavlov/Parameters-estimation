using System;
using System.IO;
using System.Security.Cryptography;

namespace lab2
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
            List<double> I = new();
            foreach (var item in File.ReadLines("txt/напряжения.txt"))
            {
                V.Add(double.Parse(item));
            }
            foreach (var item in File.ReadLines("txt/Токи.txt"))
            {
                I.Add(double.Parse(item));
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
            MFE mfe = new MFE("txt/Grid.txt", 50, 0.1, 0.01);
            mfe.BoundaryLeft = (x) => 0;
            mfe.BoundaryRight = (x) => 0;
            mfe.BoundaryTop = (x) => 0;
            mfe.BoundaryBot = (x) => 0;
            mfe.BoundaryTypeLeft = 2;
            mfe.BoundaryTypeRight = 1;
            mfe.BoundaryTypeTop = 2;
            mfe.BoundaryTypeBot = 1;
            mfe.Solve();
            MFEFunction F = new MFEFunction(A, B, M, N, I, mfe, M.Count, V);//истинные значения h=50 sigma1 =0.1 sigma2=0.01
            List<double> Params = new() { 50, 0.01 };
            var asd1 = F.CalcVi(Params, 0);
            var asd2 = F.CalcVi(Params, 1);
            var asd3 = F.CalcVi(Params, 2);
            GaussNewton solver = new GaussNewton(0.1, Params, F);
            solver.Solve(100, 1e-10);


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
            rn = text.Length-1;
            for (int i = 0; i < text.Length; i++)
            {
                rgrid.Add(double.Parse(text[i]));
            }
            zn = rn;
            for (int i = 0; i <= zn; i++)
            {
                zgrid.Add(-rgrid[zn - i]);
            }
            mat = new Matrix((rn + 1) * (zn + 1));
            Console.WriteLine($"{rgrid[1] - rgrid[0]} {rgrid[zgrid.Count - 1] - rgrid[zgrid.Count - 2]}");
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
                    double cursigma = (zgrid[z + 1] + zgrid[z]) / 2 > -h ? sigma1 : sigma2;
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
                }
            }
            mat.b[(rn + 1) * (zn)] += 1 / (2.0 * Math.PI);
        }
        private void AddBoundary()
        {

            if (BoundaryTypeBot == 1)
            {
                for (int r = 0; r <= rn; r++)
                {
                    var value = BoundaryBot(rgrid[r]);
                    mat.di[r] = 1;
                    mat.b[r] = value;
                    for (int k = mat.ia[r]; k < mat.ia[r + 1]; k++)
                    {
                        mat.b[mat.ja[k]] -= mat.al[k] * value;
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
                            mat.b[i] -= mat.al[k] * value;
                            mat.al[k] = 0;
                        }
                    }
                }
            }
            else
            {
                //2 краевые низ
            }
            if (BoundaryTypeRight == 1)
            {
                for (int z = 0; z <= zn; z++)
                {
                    var value = BoundaryRight(zgrid[z]);
                    int index = rn + (rn + 1) * z;
                    mat.di[index] = 1;
                    mat.b[index] = value;
                    for (int k = mat.ia[index]; k < mat.ia[index + 1]; k++)
                    {
                        mat.b[mat.ja[k]] -= mat.al[k] * value;
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
                            mat.b[i] -= mat.al[k] * value;
                            mat.al[k] = 0;
                        }
                    }
                }
            }
            else
            {
                //2 краевые право
            }
            if (BoundaryTypeTop == 1)
            {
                for (int r = 0; r <= rn; r++)
                {
                    var value = BoundaryTop(rgrid[r]);
                    int index = zn * (rn + 1) + r;
                    mat.di[index] = 1;
                    mat.b[index] = value;
                    for (int k = mat.ia[index]; k < mat.ia[index + 1]; k++)
                    {
                        mat.b[mat.ja[k]] -= mat.al[k] * value;
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
                            mat.b[i] -= mat.al[k] * value;
                            mat.al[k] = 0;
                        }
                    }
                }
            }
            else
            {
                //2 краевые верх
            }
            if (BoundaryTypeLeft == 1)
            {
                for (int z = 0; z <= zn; z++)
                {
                    var value = BoundaryLeft(zgrid[z]);
                    int index = (rn + 1) * z;
                    mat.di[index] = 1;
                    mat.b[index] = value;
                    for (int k = mat.ia[index]; k < mat.ia[index + 1]; k++)
                    {
                        mat.b[mat.ja[k]] -= mat.al[k] * value;
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
                            mat.b[i] -= mat.al[k] * value;
                            mat.al[k] = 0;
                        }
                    }
                }
            }
            else
            {
                //2 краевые лево
            }
        }
        private int rn;
        private int zn;
        private static class Matrices
        {
            public static double[][][] GMatr = new double[2][][]
            {
                new double[4][]
                {
                    new double[4]{5.0/12,-1.0/12,-1.0/6,-1.0/6},//x y  x
                    new double[4]{-1.0/12,1.0/4,-1.0/6,0 },//1-x y x
                    new double[4]{-1.0/6,-1.0/6,5.0/12,-1.0/12},//1-y x x
                    new double[4]{ -1.0 / 6,0,-1.0/12,1.0/4 }//1-y 1-x x
                },
                new double[4][]
                {
                    new double[4]{1.0/4,-1.0/12,0,-1.0/6},//x y 1-x
                    new double[4]{-1.0/12,5.0/12,-1.0/6,-1.0/6},//1-x y 1-x
                    new double[4]{0,-1.0/6,1.0/4,-1.0/12},
                    new double[4]{-1.0/6,-1.0/6,-1.0/12,5.0/12}
                }
            };


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
            mat.di_LU.Add(0);
            mat.b.Add(0);
            for (int i = 1; i < dim; i++)
            {
                mat.di.Add(0);
                mat.di_LU.Add(0);
                mat.b.Add(0);
                int count = 0;
                foreach (var item in list[i])
                {
                    if (item < i)
                    {
                        mat.ja.Add(item);
                        mat.al.Add(0);
                        mat.al_LU.Add(0);
                        mat.au_LU.Add(0);
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
            //Console.WriteLine($"r={r} V={q[dic[0]] * rloc * zloc + q[dic[1]] * (1 - rloc) * zloc + q[dic[2]] * rloc * (1 - zloc) + q[dic[3]] * (1 - rloc) * (1 - zloc)}");
            return q[dic[0]] * rloc * zloc + q[dic[1]] * (1 - rloc) * zloc + q[dic[2]] * rloc * (1 - zloc) + q[dic[3]] * (1 - rloc) * (1 - zloc);
        }
        private double h;
        public double H
        {
            get => h;
            set
            {
                if (h != value)
                    IsSolved = false;
                h = value;
            }
        }
        private Matrix mat;
        private double sigma1;
        public int BoundaryTypeTop = 1;
        public int BoundaryTypeBot = 1;
        public int BoundaryTypeLeft = 1;
        public int BoundaryTypeRight = 1;
        public Func<double, double> BoundaryTop = (x) => 0;
        public Func<double, double> BoundaryBot = (x) => 0;
        public Func<double, double> BoundaryLeft = (x) => 0;
        public Func<double, double> BoundaryRight = (x) => 0;
        public double Sigma1
        {
            get => sigma1;
            set
            {
                if (sigma1 != value)
                    IsSolved = false;
                sigma1 = value;
            }
        }
        private double sigma2;
        public double Sigma2
        {
            get => sigma2;
            set
            {
                if (sigma2 != value)
                    IsSolved = false;
                sigma2 = value;
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
        public List<double> di_LU = new();
        public List<double> au_LU = new();
        public List<double> al_LU = new();
        public Matrix(int n)
        {
            this.n = n;
        }
        public void Clear()
        {
            for (int i = 0; i < n; i++)
            {
                di[i] = 0;
                di_LU[i] = 0;
                b[i] = 0;
            }
            int m = al.Count;
            for (int i = 0; i < m; i++)
            {
                al[i] = 0;
                al_LU[i] = 0;
                au_LU[i] = 0;
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
        public void LU()
        {

            for (int i = 0; i < n; i++)
            {
                di_LU[i] = di[i];
            }
            var q = al.Count;
            for (int i = 0; i < q; i++)
            {
                al_LU[i] = al[i];
                au_LU[i] = al[i];
            }
            for (int i = 0; i < n; i++)
            {
                double sumdi = 0.0;

                int i0 = ia[i];
                int i1 = ia[i + 1];


                for (int k = i0; k < i1; k++)
                {
                    int j = ja[k];
                    int j0 = ia[j];

                    int j1 = ia[j + 1];


                    int ik = i0;
                    int kj = j0;

                    double suml = 0.0;
                    double sumu = 0.0;

                    while (ik < k)
                    {

                        if (ja[ik] == ja[kj])
                        {

                            suml += al_LU[ik] * au_LU[kj];
                            sumu += au_LU[ik] * al_LU[kj];
                            ik++;
                            kj++;
                        }

                        else
                        {
                            if (ja[ik] > ja[kj])
                            {
                                kj++;
                            }
                            else
                            {
                                ik++;
                            }
                        }
                    }

                    al_LU[k] = al_LU[k] - suml;
                    au_LU[k] = (au_LU[k] - sumu) / di_LU[j];
                    sumdi += al_LU[k] * au_LU[k];
                }

                di_LU[i] = di_LU[i] - sumdi;
            }
        }
        List<double> LUDirect(List<double> rpart)
        {
            List<double> res = new();
            for (int i = 0; i < n; i++)
            {
                res.Add(rpart[i]);
            }

            for (int i = 0; i < n; i++)
            {
                double sum = 0.0;
                for (int j = ia[i]; j < ia[i + 1]; j++)
                    sum += al_LU[j] * res[ja[j]];
                res[i] -= sum;
                res[i] /= di_LU[i];
            }
            return res;
        }
        List<double> LUReverse(List<double> rpart)
        {
            List<double> res = new();
            for (int i = 0; i < n; i++)
            {
                res.Add(rpart[i]);
            }

            for (int i = n - 1; i >= 0; i--)
            {
                for (int j = ia[i]; j < ia[i + 1]; j++)
                    res[ja[j]] -= au_LU[j] * res[i];
            }
            return res;
        }
        List<double> MultU(List<double> x0)
        {
            List<double> res = new();
            for (int i = 0; i < n; i++)
            {
                res.Add(x0[i]);
                int k0 = ia[i], k1 = ia[i + 1];
                for (int k = k0; k < k1; k++)
                {
                    res[i] += au_LU[k] * x0[ja[k]];
                }
            }
            return res;
        }

        public List<double> GMRES(List<double> x0, double eps, int maxiter, int depth)
        {
            //Q верхнетреугольная обратный 
            //S нижнетреугольная прямой
            int curdepth = depth;
            double bnorm = Math.Sqrt(DotProduct(b, b));
            List<double> x = MultU(x0);//меняем на x с волной

            var r = MatrixMult(x0);
            for (int i = 0; i < n; i++)
            {
                r[i] = b[i] - r[i];
            }
            r = LUDirect(r);//считаем невязку с волной

            List<List<double>> V = new();//транспонированная
            double[][] H = new double[depth + 1][];
            for (int i = 0; i < depth + 1; i++)
            {
                H[i] = new double[depth];
            }
            double rnorm = Math.Sqrt(DotProduct(r, r));
            int iter = 0;
            while (rnorm / bnorm > eps && iter < maxiter)
            {
                curdepth = depth;
                V.Add(new());
                for (int i = 0; i < n; i++)
                {
                    V[0].Add(r[i] / rnorm);// добавляем первый вектор в матрицу V
                }
                for (int i = 0; i < depth; i++)// заполняется матрица H
                {
                    List<double> w = LUReverse(V[i]);
                    w = MatrixMult(w);
                    w = LUDirect(w);
                    for (int j = 0; j <= i; j++)
                    {
                        H[j][i] = DotProduct(V[j], w);
                    }
                    for (int j = 0; j <= i; j++)
                    {
                        for (int k = 0; k < n; k++)
                        {
                            w[k] -= H[j][i] * V[j][k];
                        }
                    }

                    H[i + 1][i] = Math.Sqrt(DotProduct(w, w));
                    if (H[i + 1][i] == 0)// если новый вектор нулевой заканчиваем заполнение матриц
                    {
                        curdepth = i + 1;
                        break;
                    }
                    else
                    {
                        V.Add(new());
                        for (int k = 0; k < n; k++)
                        {
                            V[i + 1].Add(w[k] / H[i + 1][i]);//если ненулевой добавляем в матрицу V
                        }
                    }
                }
                List<double> d = new();//дальше будет минимизация d-Hz
                d.Add(rnorm);
                for (int i = 1; i < curdepth + 1; i++)
                {
                    d.Add(0);
                }
                for (int i = 0; i < curdepth; i++)//умножаем правую часть и матрицу на матрицу поворота чтобы убрать диагональ под главной
                {
                    double norm = Math.Sqrt(H[i][i] * H[i][i] + H[i + 1][i] * H[i + 1][i]);
                    double c = H[i][i] / norm;
                    double s = H[i + 1][i] / norm;
                    for (int k = i; k < curdepth; k++)
                    {
                        double ii = c * H[i][k] + s * H[i + 1][k];
                        double i1i = -s * H[i][k] + c * H[i + 1][k];
                        H[i][k] = ii;
                        H[i + 1][k] = i1i;
                    }
                    double d1 = d[i] * c + d[i + 1] * s;
                    double d2 = -s * d[i] + c * d[i + 1];
                    d[i] = d1;
                    d[i + 1] = d2;
                }

                for (int i = curdepth - 1; i >= 0; i--)//обратный гаус для матрицы H верхнетреугольной
                {
                    double summ = 0;
                    for (int j = curdepth - 1; j > i; j--)
                    {
                        summ += d[j] * H[i][j];
                    }
                    d[i] = (d[i] - summ) / H[i][i];
                }

                for (int i = 0; i < n; i++)//делаем добавку в x
                {
                    for (int j = 0; j < curdepth; j++)
                    {
                        x[i] += V[j][i] * d[j];
                    }
                }

                r = LUReverse(x);
                r = MatrixMult(r);
                for (int i = 0; i < n; i++)
                {
                    r[i] = b[i] - r[i];
                }
                r = LUDirect(r);
                rnorm = Math.Sqrt(DotProduct(r, r));
                iter++;
                V.Clear();
             }
            Console.WriteLine($"{iter} {rnorm / bnorm}");
            return LUReverse(x);
        }

        public List<double> LoS_precond(List<double> x0, double eps, int maxiter)
        {
            int k = 1;
            List<double> buf = MatrixMult(x0);
            double bnorm = 0;
            for (int i = 0; i < n; i++)
            {
                buf[i] = b[i] - buf[i];
            }
            double rnorm = Math.Sqrt(DotProduct(buf, buf));
            List<double> r = LUDirect(buf);
            bnorm = Math.Sqrt(DotProduct(b, b));
            List<double> z = LUReverse(r);
            buf = MatrixMult(z);
            List<double> p = LUDirect(buf);
            double resid = 1;
            while (resid > eps && k < maxiter)
            {
                double pp = DotProduct(p, p);
                double pr = DotProduct(p, r);
                double alpha = pr / pp;
                for (int i = 0; i < n; i++)
                {
                    x0[i] += alpha * z[i];
                    r[i] -= alpha * p[i];
                }
                rnorm = Math.Sqrt(DotProduct(r, r));
                List<double> Ur = LUReverse(r);
                buf = MatrixMult(Ur);
                buf = LUDirect(buf);
                double betta = -(DotProduct(p, buf) / pp);
                for (int i = 0; i < n; i++)
                {
                    z[i] = Ur[i] + betta * z[i];
                    p[i] = buf[i] + betta * p[i];
                }
                double test1 = 0;
                double test2 = 0;
                var asd = MatrixMult(x0);
                for (int i = 0; i < n; i++)
                {
                    test1 += (asd[i] - b[i]) * (asd[i] - b[i]);
                    test2 += b[i] * b[i];
                }
                resid = Math.Sqrt(test1 / test2);
                k++;
            }
            Console.WriteLine($"{k} {rnorm / bnorm} {resid}");
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
            LU();
            return GMRES(x0, 1e-12, 1000,30);
        }
    }
    public class MatrixFull
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
        public MatrixFull(int n)
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
    public abstract class Function
    {
        public int n;//количество точек в которых были измерения
        public List<double> V;//Значения измерений

        protected Function(int n, List<double> v)
        {
            this.n = n;
            V = v;
        }

        public abstract double Calc(List<double> Params);
        public abstract double CalcVi(List<double> Params, int i);
    }
    public class MFEFunction : Function
    {
        List<Vector> A;
        List<Vector> B;
        List<Vector> M;
        List<Vector> N;
        List<double> I;
        public MFE mfe;

        public MFEFunction(List<Vector> A, List<Vector> B, List<Vector> M, List<Vector> N, List<double> I, MFE mfe, int n, List<double> V) : base(n, V)
        {
            this.A = A;
            this.B = B;
            this.M = M;
            this.N = N;
            this.I = I;
            this.mfe = mfe;
        }

        public override double Calc(List<double> Params)
        {
            mfe.H = Params[0];
            mfe.Sigma2 = Params[1];
            double res = 0;
            int priem = M.Count;
            int istok = A.Count;
            for (int i = 0; i < priem; i++)
            {
                double curV = 0;
                for (int j = 0; j < istok; j++)
                {
                    curV += I[j] * (mfe.GetSollution((B[j] - M[i]).Norm, 0) - mfe.GetSollution((A[j] - M[i]).Norm, 0) - mfe.GetSollution((B[j] - N[i]).Norm, 0) + mfe.GetSollution((A[j] - N[i]).Norm, 0));
                }
                res += ((curV - V[i]) / V[i]) * ((curV - V[i]) / V[i]);
            }
            return res;
        }

        public override double CalcVi(List<double> Params, int i)
        {
            mfe.Sigma2 = Params[1];
            mfe.H = Params[0];
            int istok = A.Count;
            double curV = 0;
            for (int j = 0; j < istok; j++)
            {
                curV += I[j] * (mfe.GetSollution((B[j] - M[i]).Norm, 0) - mfe.GetSollution((A[j] - M[i]).Norm, 0) - mfe.GetSollution((B[j] - N[i]).Norm, 0) + mfe.GetSollution((A[j] - N[i]).Norm, 0));
            }
            return curV;
        }
    }
    public class GaussNewton
    {
        public double diffparameter;
        public List<double> Params;
        Function F;

        public GaussNewton(double diffparameter, List<double> @params, Function f)
        {
            this.diffparameter = diffparameter;
            Params = @params;
            F = f;
        }

        public List<double> Solve(int maxiter, double eps)
        {
            List<double> weights = new();
            for (int i = 0; i < F.n; i++)
            {
                weights.Add(1 / F.V[i]);
            }
            int n = Params.Count;
            MatrixFull mat = new MatrixFull(n);
            List<double> b = new();
            for (int i = 0; i < n; i++)
            {
                b.Add(0);
            }
            double penalty = F.Calc(Params);
            int k = 0;
            while (penalty > eps && k < maxiter)
            {
                double[][] diffs = new double[F.n][];
                double[] V0 = new double[F.n];
                for (int i = 0; i < F.n; i++)
                {
                    V0[i] = F.CalcVi(Params, i);
                    diffs[i] = new double[n];
                }

                for (int i = 0; i < n; i++)
                {
                    double delta = Params[i] * diffparameter;
                    Params[i] += delta;
                    for (int s = 0; s < F.n; s++)
                    {
                        diffs[s][i] = F.CalcVi(Params, i);
                        diffs[s][i] = (diffs[s][i] - V0[s]) / delta;
                    }
                    Params[i] -= delta;

                }

                for (int s = 0; s < F.n; s++)
                {
                    for (int i = 0; i < n; i++)
                    {
                        b[i] += weights[s] * weights[s] * diffs[s][i] * (F.V[s] - V0[s]);
                        for (int j = 0; j < n; j++)
                        {
                            mat.mat[i][j] += weights[s] * weights[s] * diffs[s][i] * diffs[s][j];
                        }
                    }
                }
                for (int i = 0; i < n; i++)
                {
                    mat.mat[i][i] += 1e-5;
                }
                //добавить регуляризацию
                var dparam = mat.SolveLU(b);
                mat.Clear();
                for (int i = 0; i < n; i++)
                {
                    Params[i] += dparam[i];
                    b[i] = 0;
                }
                penalty = F.Calc(Params);
                Console.WriteLine($"{penalty} {(F as MFEFunction).mfe.Sigma2 } {(F as MFEFunction).mfe.H }");
                k++;
            }
            return Params;
        }
    }
}