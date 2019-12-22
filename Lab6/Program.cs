using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static System.Console;
using static System.Math;

namespace Lab6
{
    class Program
    {
        static void Main(string[] args)
        {
            WriteLine("Введите номер метода:\n1) Метод наискорейшего спуска\n2) Метод Давидона-Флетчера-Пауэлла\n3) Метод Флетчера-Ривса ");
            int n = int.Parse(ReadLine());
            switch(n)
            {
                case 1:
                    SD_Call();
                    break;
                case 2:
                    DFP_Call();
                    break;
                case 3:
                    FR_Call();
                    break;
                default:
                    return;
                
            }
            ReadKey();
        }
        private static void FR_Call()
        {
            OutputEncoding = System.Text.Encoding.UTF8;
            ForegroundColor = ConsoleColor.Green;
            FR_Basic.Start(new double[] { 1, 2, 3 }, 3, FR_Basic.F_Main, FR_Basic.G_Main);
            ForegroundColor = ConsoleColor.White;
        }

        private static void DFP_Call()
        {
            OutputEncoding = System.Text.Encoding.UTF8;
            ForegroundColor = ConsoleColor.Green;
            WriteLine("Выберите функцию\n 1) (x1 + 10*x2)^2 + 5*(x3-x4)^2 + (x2-2*x3)^4 + 10*(x1-x4)^4\n 2)3*(x1-4)^2 + 5*(x2+3)^2 + 7*(2*x3 + 1)");
            int C = int.Parse(ReadLine());
            if (C == 1) DFP_Basic.Start(new double[] { -144, 67, 4 }, 3,
                DFP_Basic.F_Main, DFP_Basic.G_Main);
            else if (C == 2) DFP_Basic.Start(new double[] { 3, -1, 0, 1 },
                4, DFP_Basic.F_Bandy, DFP_Basic.G_Bandy);
            else
            {
                WriteLine();
                DFP_Call();
            }
            ForegroundColor = ConsoleColor.White;
        }

        private static void SD_Call()
        {

            OutputEncoding = System.Text.Encoding.UTF8;

            ForegroundColor = ConsoleColor.Green;

            WriteLine("Метод наискорейшего спуска");

            WriteLine($"Выберите функцию\n1)(x1 - 1)^2 + (x2-3)^2 + 4*(x3+5)^2\n2)(x1 - 1)^4 + (x2-3)^2 + 4*(x3+5)^4");

            Boolean para = true;

            if (int.TryParse(ReadLine(), out int C) && C > 0 && C < 3)

                para = C == 1;

            List<double> x0 = new List<double>();

            WriteLine("Введите начальные точки:");

            for (int i = 0; i < 3; i++) if (double.TryParse(ReadLine(), out double V)) x0.Add(V);

                else WriteLine("Значение не принято, введите его ещё раз");

            Functions Parametres = new Functions(para);

            List<double> res = SteepestDescentMethod.GradientDescent(x0, Parametres);
            if (para == true)
            {
                WriteLine($"Минимум функции (x1 - 1)^2 + (x2-3)^2 + 4*(x3+5)^2 = {Parametres.Function(res).ToString()}");
            }
            else
            {
                WriteLine($"Минимум функции (x1 - 1)^4 + (x2-3)^2 + 4*(x3+5)^4 = {Parametres.Function(res).ToString()}");
            }
            

            WriteLine("Точки минимума:");

            for (int i = 0; i < res.Count; i++) Write($"x{(Char)(8321 + i)} = {res[i]};\t");

            GC.Collect();

            ForegroundColor = ConsoleColor.White;
        }



    }

    /// <summary>
    /// Параметры задач - функции, частные производные и градиент
    /// </summary>
    public class Functions
    {
        #region Поля и конструктор
        public Func<List<double>, List<double>> Gradient;
        public Func<List<double>, double> Function;
        public Func<double, int, double> PartialFunction;

        /// <summary>
        /// True - задание 1<para/>
        /// False - задание 2
        /// </summary>
        /// <param name="N"></param>
        public Functions(bool N)
        {
            if (N)
            {
                Gradient = GradientF1;
                Function = Function1;
                PartialFunction = PartialFunction1;
            }
            else
            {
                Gradient = GradientF2;
                Function = Function2;
                PartialFunction = PartialFunction2;
            }
        }
        #endregion

        #region Задание 1
        private static List<double> GradientF1(List<double> X)
        {
            List<double> tmp = new List<double>();
            // Частные производные
            tmp.Add(2 * (X[0] - 1));
            tmp.Add(2 * (X[1] - 3));
            tmp.Add(8 + (X[2] + 5));
            return tmp;
        }

        private static double Function1(List<double> X) =>
            Math.Pow(X[0] - 1, 2) + Math.Pow(X[1] - 3, 2) + 4 * Math.Pow(X[2] + 5, 2);

        private static double PartialFunction1(double x, int i)
        {
            double z;
            switch (i)
            {
                case 0:
                    z = Math.Pow(x - 1, 4);
                    break;
                case 1:
                    z = Math.Pow(x - 3, 2);
                    break;
                case 2:
                    z = Math.Pow(x + 5, 4) * 4;
                    break;
                default:
                    z = 0;
                    break;
            }
            return z;
        }
        #endregion

        #region Задание 2
        private static List<double> GradientF2(List<double> X)
        {
            List<double> tmp = new List<double>();
            // Частные производные
            tmp.Add(4 * (X[0] - 1));
            tmp.Add(2 * (X[1] - 3));
            tmp.Add(16 + (X[2] + 5));
            return tmp;
        }

        private static double Function2(List<double> X) =>
            Math.Pow(X[0] - 1, 4) + Math.Pow(X[1] - 3, 2) + 4 * Math.Pow(X[2] + 5, 4);

        private static double PartialFunction2(double x, int i)
        {
            double z;
            switch (i)
            {
                case 0:
                    z = Math.Pow(x - 1, 4);
                    break;
                case 1:
                    z = Math.Pow(x - 3, 2);
                    break;
                case 2:
                    z = Math.Pow(x + 5, 4) * 4;
                    break;
                default:
                    z = 0;
                    break;
            }
            return z;

        }
        #endregion

    }

    public partial class SteepestDescentMethod
    {
        private static List<double> GoldenSection_v1(List<double> x, double a, double b, int n)
        {
            // Параметр n здесь не исполюзуется, но удалять
            // я его не стал для сохранения однообразия функций
            List<double> Result = new List<double>();
            for (int i = 0; i < x.Count; i++)
                Result.Add(GS_v1_Iter(a, b, i));
            return Result;
        }

        private static double GS_v1_Iter(double A, double B, int i)
        {
            double T1, T2, X;
            T1 = 0.3819660113;
            T2 = 1 - T1;
            double X0, X1, X2, X3;
            X0 = A; X1 = A + T1 * (B - A);
            X2 = A + T2 * (B - A); X3 = B;
            X = X1;
            double F1, F2, I;
            F1 = Parametres.PartialFunction(X, i);
            X = X2; F2 = Parametres.PartialFunction(X, i);
            do
            {
                if (F1 < F2)
                {
                    I = X2 - X0; X3 = X2; X2 = X1; X1 = X0 + T1 * I;
                    F2 = F1; X = X1; F1 = Parametres.PartialFunction(X, i);
                }
                else
                {
                    I = X3 - X1; X0 = X1; X1 = X2; X2 = X0 + T2 * I;
                    F1 = F2; X = X2; F2 = Parametres.PartialFunction(X, i);
                }
            } while (I > 0.001);
            return X1;
        }

        #region Неудачная попытка избежать цикла
        private static List<double> GoldenSection_v2(List<double> x, double a, double b, int n)
        {
            List<double> tmp = x;
            int i, j;
            double s1;
            double s2;
            double u1;
            double u2;
            double fu1;
            double fu2;
            List<double> GF = Parametres.Gradient(x);

            s1 = (3 - Math.Sqrt(5)) / 2;
            s2 = (Math.Sqrt(5) - 1) / 2;
            u1 = a + s1 * (b - a);
            u2 = a + s2 * (b - a);

            for (j = 0; j < x.Count; j++)
                tmp[j] = x[j] + u1 * GF[j];
            fu1 = Parametres.Function(tmp);

            for (j = 0; j < x.Count; j++)
                tmp[j] = x[j] + u2 * GF[j];
            fu2 = Parametres.Function(tmp);

            for (i = 1; i <= n; i++)
            {
                if (fu1 <= fu2)
                {
                    b = u2;
                    u2 = u1;
                    fu2 = fu1;
                    u1 = a + s1 * (b - a);
                    for (j = 0; j < x.Count; j++)
                        tmp[j] = x[j] + u1 * GF[j];
                    fu1 = Parametres.Function(tmp);
                }
                else
                {
                    a = u1;
                    u1 = u2;
                    fu1 = fu2;
                    u2 = a + s2 * (b - a);
                    for (j = 0; j < x.Count; j++)
                        tmp[j] = x[j] + u2 * GF[j];
                    fu2 = Parametres.Function(tmp);
                }
            }

            for (j = 0; j < x.Count; j++)
                tmp[j] = x[j] + u1 * GF[j];

            fu1 = Parametres.Function(tmp);

            for (j = 0; j < x.Count; j++)
                tmp[j] = x[j] + u2 * GF[j];

            fu2 = Parametres.Function(tmp);

            if (fu1 < fu2) for (j = 0; j < x.Count; j++)
                    tmp[j] = x[j] + u1 * GF[j];

            else for (j = 0; j < x.Count; j++)
                    tmp[j] = x[j] + u2 * GF[j];

            return tmp;
        }
        #endregion

    }


    public static partial class SteepestDescentMethod
    {
        static double MaxIterCount = 100000;

        static double EPS = Math.Pow(10, -5);

        static Func<List<double>, double, double, int, List<double>> GetLinearMin;

        static Functions Parametres;

        /// <summary>
        /// Метод наискорейшего спуска
        /// </summary>
        /// <param name="x0">Переменные</param>
        /// <returns></returns>
        public static List<double> GradientDescent(List<double> x0, Functions _Parametres)
        {
            Parametres = _Parametres;
            List<double> old = new List<double>(), cur_x = x0, gr;
            //Console.ForegroundColor = ConsoleColor.Green;

            // Выбор метода линейного поиска
            Console.WriteLine("Выберите метод линейного поиска" +
                              "\n1) Метод квадратической интерполяции" +
                              "\n2) Метод Фибоначчи" +
                              "\n3) Метод \"Золотого Сечения\"");
            int Choise = int.Parse(Console.ReadLine());
            switch (Choise)
            {
                case 1:
                    GetLinearMin = QuadraticInterpolation;
                    break;
                case 2:
                    GetLinearMin = Fibonacci;
                    break;
                case 3:
                    GetLinearMin = GoldenSection_v1;
                    break;
            }

            for (int Iterations = 1; Iterations <= MaxIterCount; Iterations++)
            {
                for (int y = 0; y < x0.Count; y++) old.Add(cur_x[y]);
                gr = Parametres.Gradient(cur_x);

                #region Линейный поиск
                cur_x = GetLinearMin(cur_x, -10, 10, 100);
                #endregion

                //условие останова
                double s = Math.Abs(Parametres.Function(cur_x) - Parametres.Function(old));
                if (s < EPS) return cur_x;
            }
            return cur_x;
        }


    }




    public partial class SteepestDescentMethod
    {
        private static List<double> QuadraticInterpolation(List<double> x, double a, double b, int n)
        {
            List<double> Result = new List<double>();
            for (int i = 0; i < x.Count; i++)
                Result.Add(GetMinimum((a + b) / 2, 0.5, Parametres.PartialFunction, i));
            return Result;
        }

        private static Double GetMinimum(Double A, Double H, Func<Double, int, Double> func, int ii)
        {
            double _X;
            double E = 0.000001;
            double[] X = new double[4];
            double[] F = new double[4];
            X[0] = A;
            _X = X[0];
            F[0] = func(_X, ii);
            X[1] = A + H;
            _X = X[1];
            F[1] = func(_X, ii);
            if (F[0] < F[1])
            {
                X[2] = A - H; _X = X[2]; F[2] = func(_X, ii);
            }
            else
            {
                X[2] = A + 2 * H; _X = X[2]; F[2] = func(_X, ii);
            }
            //вычисление первого аппроксимируещего минимума 
            double DN, NM;
            DN = (X[1] - X[2]) * F[0];
            DN = DN + (X[2] - X[0]) * F[1] + (X[0] - X[1]) * F[2];
            NM = (X[1] * X[1] - X[2] * X[2]) * F[0];
            NM = NM + (X[2] * X[2] - X[0] * X[0]) * F[1];
            NM = NM + (X[0] * X[0] - X[1] * X[1]) * F[2];
            X[3] = NM / (2 * DN);
            _X = X[3]; F[3] = func(_X, ii);
            while (true)
            {
                double _F;
                //упорядочить значение функции
                for (int J = 0; J < 4; J++)
                {
                    for (int K = J + 1; K < 4; K++)
                    {
                        if (F[J] <= F[K])
                            continue;
                        _X = X[J]; X[J] = X[K]; X[K] = _X;
                        _F = F[J]; F[J] = F[K]; F[K] = _F;
                        /*поменять местами F[J] и F[K],а также X[J] и X[K],
                         * если они не упорядочены*/
                    }
                }
                double S1, S2, S3;
                //закончить,если получена заданная точность
                if (Math.Abs(X[0] - X[1]) >= E)
                {
                    //запомнить три лучшие точки
                    S1 = Math.Sign(X[0] - X[1]);
                    S2 = Math.Sign(X[2] - X[0]);
                    S3 = Math.Sign(X[3] - X[0]);
                    if (S1 == S2 && S1 == -S3)
                    {
                        X[2] = X[3]; F[2] = F[3];
                    }
                    //вторая интерполяция 
                    DN = (X[1] - X[2]) * F[0] + (X[2] - X[0]) * F[1] + (X[0] - X[1]) * F[2];
                    _F = (F[0] - F[1]) / (2 * DN);
                    _F = _F * (X[1] - X[2]) * (X[2] - X[0]);
                    X[3] = (X[0] + X[1]) / 2 + _F;
                    _X = X[3]; F[3] = func(_X, ii);
                    //повторить вторую интерполяцию 
                }
                else break;
            }
            return X[0];
        }
    }


    




    public partial class SteepestDescentMethod
    {
        private static List<double> Fibonacci(List<double> x, double a, double b, int n)
        {
            List<double> Result = new List<double>();
            for (int i = 0; i < x.Count; i++)
                Result.Add(FibonachiRun(a, b, n, i));
            return Result;
        }

        private static double FibonachiRun(double a, double b, int n, int i)
        {
            int N = 0;
            int[] Fib = new int[2 * n + 1];
            Fib[0] = Fib[1] = 1;

            #region Вычисление чисел Фибоначчи
            for (int j = 2; j <= 2 * n; j++)
                Fib[j] = Fib[j - 1] + Fib[j - 2];
            #endregion

            double X1, X2, F2, F4;
            int G = 0;
            double FN1 = 1, FN2 = 1, FN,
                F = (b - a) / EPS;

            while (FN1 < F)
            {
                FN = FN1 + FN2;
                FN1 = FN2;
                FN2 = FN;
                N++;
            }
            bool bix;
            int ix = N & 1;
            if (ix == 1)
                bix = true;
            else
                bix = false;
            X1 = a + (double)Fib[N - 2] / Fib[N] * (b - a) - (bix ? -1 : 1) * EPS / Fib[N];
            X2 = a + (double)Fib[N - 1] / Fib[N] * (b - a) + (bix ? -1 : 1) * EPS / Fib[N];
            F2 = Parametres.PartialFunction(X1, i);
            F4 = Parametres.PartialFunction(X2, i);
            while (Math.Abs(b - a) > EPS)
            {
                if (F2 >= F4)
                {
                    ix = (N - G) & 1;
                    if (ix == 1)
                        bix = true;
                    else
                        bix = false;
                    a = X1;
                    X1 = X2;
                    F2 = F4;
                    X2 = a + (double)Fib[N - G - 1] / Fib[N - G] * (b - a) + (bix ? -1 : 1) * EPS / Fib[N - G];
                    F4 = Parametres.PartialFunction(X2, i);
                }
                else
                {
                    ix = (N - G) & 1;
                    if (ix == 1)
                        bix = true;
                    else
                        bix = false;
                    b = X2;
                    X2 = X1;
                    F4 = F2;
                    X1 = a + (double)Fib[N - G - 2] / Fib[N - G] * (b - a) - (bix ? -1 : 1) * EPS / Fib[N - G];
                    F2 = Parametres.PartialFunction(X1, i);
                }
                G++;
            }
            return (a + b) / 2;
        }
    }

   
    public static class DFP_Basic
    {
        static double[] X;
        static double Z;
        static double[] G;
        static double G0;
        public static void Start(double[] x, int N, Action F, Action G6)
        {
            double FP, GP, QX, G1, HH, BB, FQ, G2, GQ, ZZ, WW, W, DD, FR, GR, G3, KK, DK, WK;
            X = x; double CC = 0;//CC-колво итераций
            double[] P = new double[N], V = new double[N], Y = new double[N], M = new double[N], U = new double[N], Q = new double[N], D = new double[N]; G = new double[N];
            double[,] H = new double[N, N];
            for (int i = 0; i < N; i++)
            {
                H[i, i] = 1;
            }
            do
            {
                for (int i = 0; i < N; i++)
                {
                    P[i] = X[i]; Y[i] = X[i]; WriteLine($"X {i + 1} {X[i]}");
                }
                WriteLine();
                WriteLine();
                F(); WriteLine($"[{CC}] Z = {Z}"); //420
                FP = Z; G6(); G1 = G0;
                //Градиент запомнить в u и выбрать начальное значение d
                for (int i = 0; i < N; i++)
                {
                    U[i] = G[i]; D[i] = 0;
                    for (int j = 0; j < N; j++)
                    {
                        D[i] = D[i] - H[i, j] * G[j];
                    }
                }
                while (true)
                {
                    GP = 0;
                    for (int i = 0; i < N; i++) GP = GP + G[i] * D[i];//610
                    if (GP < 0) break;
                    // найти начальный шаг если необходимо
                    // измеить направление спуска на противоположное
                    QX = Abs(2 * FP / GP); if (QX > 1) QX = 1;
                    for (int i = 0; i < N; i++)
                    {
                        X[i] = P[i] - QX * D[i]; P[i] = X[i];
                    }
                    F(); FP = Z; WriteLine("Нестабильность?");
                    G6(); G1 = G0;
                }
                QX = Abs(2 * FP / GP); if (QX > 1) QX = 1;//680
                HH = QX;
                while (true)
                {
                    //найти следующую точку Q
                    BB = HH;
                    for (int i = 0; i < N; i++)
                    {
                        Q[i] = P[i] + BB * D[i]; X[i] = Q[i];
                    }
                    F(); FQ = Z;
                    G6(); G2 = G0;
                    GQ = 0;
                    for (int i = 0; i < N; i++) GQ = GQ + G[i] * D[i];
                    if (GQ > 0 || FQ > FP) break; //goto 830
                    HH = 2 * HH;
                }
                while (true)
                {
                    ZZ = 3 * (FP - FQ) / HH; ZZ = ZZ + GP + GQ;//830
                    WW = ZZ * ZZ - GP * GQ; if (WW < 0) WW = 0;
                    W = Sqrt(WW);
                    DD = HH * (1 - (GQ + W - ZZ) / (GQ - GP + 2 * W));
                    for (int i = 0; i < N; i++) X[i] = P[i] + DD * D[i];//870
                    F(); FR = Z;
                    G6(); G3 = G0;
                    //найти градиент в новой точке
                    GR = 0;
                    for (int i = 0; i < N; i++) GR = GR + G[i] * D[i];//910
                    if ((Z <= FP && Z <= FQ) || GR > 0) break; //goto 990
                    HH = HH - DD;
                    for (int i = 0; i < N; i++) P[i] = X[i];//970
                    FP = Z; GP = GR; G1 = G0;
                }
                HH = DD;//990
                if (!(Z <= FP && Z <= FQ)) // if (Z<=FP && Z<=FQ) goto 1100
                {
                    HH = DD;
                    //обновить матрцу H
                }
                KK = 0; WK = 0; DK = 0;//1100
                for (int i = 0; i < N; i++)
                {
                    U[i] = G[i] - U[i]; V[i] = X[i] - Y[i];
                }
                for (int i = 0; i < N; i++)
                {
                    M[i] = 0;
                    for (int j = 0; j < N; j++)
                    {
                        M[i] = M[i] + H[i, j] * U[j];
                    }
                    KK = KK + M[i] * U[i]; WK = WK + V[i] * U[i];
                    DK = DK + V[i] * V[i];
                }
                //1205
                if (!(KK == 0 || WK == 0))
                {
                    for (int i = 0; i < N; i++)
                        for (int j = 0; j < N; j++)
                            H[i, j] = H[i, j] - M[i] * M[j] / KK + V[i] * V[j] / WK;
                }
                CC++;
                //проверка критерия завершения
            } while (!(Sqrt(DK) < 0.00005 || G3 < 0.00001));
            WriteLine("Минимум найден");
            WriteLine($"Количество итераций={CC} Значение минимума= {Z}");
            for (int i = 0; i < N; i++)
                WriteLine($"X {i + 1} = {X[i]}  ");
            WriteLine();//1350 end
        }

        #region Bandy Example
        public static void F_Bandy() //5000
        {
            Z = Pow(X[0] + 10 * X[1], 2) + 5 * Pow(X[2] - X[3], 2) + Pow(X[1] - 2 * X[2], 4) + 10 * Pow(X[0] - X[3], 4);
        }

        public static void G_Bandy() //6000
        {
            G[0] = 2 * (X[0] + 10 * X[1]) + 40 * Pow(X[0] - X[3], 3);
            G[1] = 20 * (X[0] + 10 * X[1]) + 4 * Pow(X[1] - 2 * X[2], 3);
            G[2] = 10 * (X[2] - X[3]) - 8 * Pow(X[1] - 2 * X[2], 3);
            G[3] = -10 * (X[2] - X[3]) - 40 * Pow(X[0] - X[3], 3);
            for (int i = 0; i < 4; i++)
            {
                G0 += G[i] * G[i];
            }
            G0 = Sqrt(G0);
            //later
        }
        #endregion

        #region Stupid Example
        public static void F_Stupid()
        {
            Z = Pow(X[0] - 1, 2) + 2.25;
        }

        public static void G_Stupid()
        {
            G[0] = 2 * X[0] - 2;
            G0 = G[0];
        }
        #endregion

        #region Main Example
        public static void F_Main()
        {
            Z = 3 * Pow(X[0] - 4, 2) + 5 * Pow(X[1] + 3, 2) + 7 * Pow(2 * X[2] + 1, 2);
        }

        public static void G_Main()
        {
            G[0] = 6 * X[0] - 24;
            G[1] = 10 * X[1] + 30;
            G[2] = 56 * X[2] + 28;
            for (int i = 0; i < 3; i++)
            {
                G0 += G[i] * G[i];
            }
            G0 = Sqrt(G0);
        }
        #endregion
    }

    public static class FR_Basic
    {
        static double[] G;
        static double[] X;
        static double Z;
        static double G0;
        public static void Start(double[] x, int N, Action F, Action G6)
        {
            bool bl;
            X = x;
            //
            int sv = 1; int tv = 0;
            double FP, GK, G1, G2, G3, K, GP, QX, HH, BB, ZZ, WW, AK, FQ, GQ, W, SV = 0, DV = 0, DD, GR, FR;
            double[] Y = new double[N], P = new double[N], Q = new double[N], D = new double[N];
            G = new double[N];
            //
            void Restart()
            {
                for (int i = 0; i < N; i++) { P[i] = X[i]; WriteLine($"X{i} = {X[i]}"); }//550
                F(); FP = Z; WriteLine($"Z = {Z}");
                G6(); G1 = G0; GK = G0;
                //В качестве первого направления взять
                // Направление наискорейшего спуска
                for (int i = 0; i < N; i++) D[i] = -G[i];
                // K - Счетчик итераций
                K = 1;
            } //[550;600)



            WriteLine("Минимизация методом Флетчера-Ривса");
            // одномерный поиск производиться кубической интерполяцией
            // Промежуточный вывод
            WriteLine("Текущие значения");
            Restart();
            do
            {
                while (true)
                {
                    GP = 0; // 600
                    for (int i = 0; i < N; i++) GP = GP + G[i] * D[i];
                    if (GP <= 0) break;
                    // Определить начальный шаг и если необходимо сменить направление спуска...
                    QX = Abs(2 * FP / GP); if (QX > 1) QX = 1;
                    for (int i = 0; i < N; i++)
                    { X[i] = P[i] - QX * D[i]; P[i] = X[i]; }
                    F(); FP = Z; WriteLine("Нестабильность!");
                    G6(); G1 = G0;
                }
                QX = Abs(2 * FP / GP); if (QX > 1) QX = 1;//680
                HH = QX;
                // Найти след точку
                while (true)
                {
                    BB = HH; // 710
                    for (int i = 0; i < N; i++)
                    {
                        Q[i] = P[i] + BB * D[i]; X[i] = Q[i];
                    }
                    F(); FQ = Z;
                    G6(); G2 = G0;
                    GQ = 0;
                    for (int i = 0; i < N; i++) GQ = GQ + G[i] * D[i];//780 - 800
                    if (GQ > 0 || FQ > FP) break;
                    //Выполнить интерполяцию
                    //удвоить шаг и перейти к точке Q
                    // или удвоить шаг
                    HH = 2 * HH;
                    for (int i = 0; i < N; i++) P[i] = Q[i];//830
                    FP = FQ; GP = GQ; G1 = G2;
                }
                while (true)
                {
                    ZZ = 3 * (FP - FQ) / HH; ZZ = ZZ + GP + GQ;// 860
                    WW = ZZ * ZZ - GP * GQ; if (WW < 0) WW = 0;
                    W = Sqrt(WW);//880
                    DD = HH * (1 - (GQ + W - ZZ) / (GQ - GP + 2 * W));
                    for (int i = 0; i < N; i++) X[i] = P[i] + DD * D[i];
                    F(); FR = Z;
                    G6(); G3 = G0;
                    //найти градиент в новой точке
                    GR = 0;
                    for (int i = 0; i < N; i++) GR = GR + G[i] * D[i];//940
                    if (Z <= FP && Z <= FQ) break;//goto 1100
                    if (GR > 0)
                    {
                        HH = DD;
                        for (int i = 0; i < N; i++) Q[i] = X[i];
                        FQ = Z; GP = GR;
                    }
                    else
                    {
                        HH = HH - DD;
                        for (int i = 0; i < N; i++) P[i] = X[i];
                        FP = Z; GP = GR;
                    }
                }
                // Проверка критерия завершения
                bl = !(G3 < .00001);
                if (bl)
                {
                    if (K == N)
                    {
                        WriteLine("Рестарт"); SV++; DV++;
                        WriteLine($"   Итерация:{SV}  Поиск:{DV}");
                        WriteLine();
                        Restart();
                    }
                    else
                    {
                        // счетичк
                        K++;
                        //Найти сопряженное направление
                        AK = G3 * G3 / (GK * GK);
                        for (int i = 0; i < N; i++) { D[i] = -G[i] + AK * D[i]; P[i] = X[i]; }
                        DV++;
                        WriteLine($"Новое направление поиска {DV}");
                        FP = Z; G1 = G0; GK = G0;
                        for (int i = 0; i < N; i++) WriteLine($"X{i} = {X[i]}"); WriteLine($"Z = {Z}");
                    }
                }
                WriteLine($"G3 = {G3}  K = {K}");
                ReadKey();
            }
            while (bl);       //1100 // return to 600

            WriteLine("Минимум найден");//1300
            for (int i = 0; i < N; i++) WriteLine($"X{i} = {X[i]}");
            WriteLine($"Минимум функции равен {Z}");
        }

        #region Main Example
        public static void F_Main()
        {
            Z = 3 * Pow(X[0] - 4, 2) + 5 * Pow(X[1] + 3, 2) + 7 * Pow(2 * X[2] + 1, 2);
        }

        public static void G_Main()
        {
            G0 = 0;
            G[0] = 6 * X[0] - 24;
            G[1] = 10 * X[1] + 30;
            G[2] = 56 * X[2] + 28;
            for (int i = 0; i < 3; i++)
            {
                G0 += G[i] * G[i];
            }
            G0 = Sqrt(G0);
        }
        #endregion
    }





}















