using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Lab2v2
{
    internal class MatrixFull
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
}
