#include "utils.h"
#include <vector>
#include <string>

SystemState gomori_solve(const Matrix& A, const Matrix& B,
    const std::vector<long double>& c) {
    SystemState ss = simplex_init(A, B, c);
    int max_iterations = 100;
    int iter = 0;
    while (true) {
        //std::cout << "#" << iter << std::endl;
        std::cout << "#" << std::left << std::setfill('#') << std::setw(c.size() * 9 + 32) << " " + std::to_string(iter) + " " << std::endl << std::setfill(' ');
        ss = simplex_step(ss);

        std::cout << "Current solution:" << std::endl;
        for (int i = 0; i < ss.opt_vars; i++) {
            long double err = abs(ss.x_curr[i] - round(ss.x_curr[i]));
            std::cout << ss.x_curr[i];
            std::cout << "(" << err << ");";
        }
        std::cout << std::endl << std::endl;

        if (ss.result != FINITE_SOL) {
            return ss;
        }
        long double tolerance = 1e-5;
        bool is_int = true;
        for (int i = 0; is_int && i < ss.opt_vars; i++) {
            long double err = abs(ss.x_curr[i] - round(ss.x_curr[i]));
            is_int = err < tolerance;
        }
        if (is_int)
            return ss;
        if (++iter >= max_iterations) {
            //std::cout << "ITERATION LIMIT!" << std::endl;
            ss.result = ITER_LIMIT;
            return ss;
        }

        int m_prev = ss.A.rows();
        int n_1_prev = ss.A.cols();
        int n_prev = n_1_prev - m_prev;
        // std::cout << "tut1" << std::endl;

        // добавляем переменную и ограничение
        ss = push_variable(ss);

        // считаем дельта, чтобы выбрать новый базисный вектор
        int n_1_curr = ss.A.cols();
        std::vector<long double> z(n_1_curr, 0.);
        std::vector<long double> Delta(n_1_curr, 0.);
        calc_table(ss.A, ss.c_ext, z, Delta);

        // std::cout << "tut3" << std::endl;
        int new_bas_idx = -1;
        //print_table(ss.A, ss.B, ss.c_ext, z, Delta, z, ss.order);

        for (int i = m_prev; i < n_1_prev; i++) {
            //std::cout << Delta[i] << " ";
            if (Delta[i] <= 0)
                continue;

            if (new_bas_idx < 0 || Delta[i] < Delta[new_bas_idx]) {
                new_bas_idx = i;
            }
        }
        //std::cout << std::endl;

        // std::cout << "tut4" << std::endl;
        // assert new_bas_idx >= m_prev
        auto A_can = ss.A;
        auto B_can = ss.B;
        auto c_ext = ss.c_ext;
        auto order = ss.order;
        // сдвигаем базисный вектор к остальным базисным
        swap_cols(A_can, c_ext, new_bas_idx, m_prev);
        // запоминаем какой переменной соответсвует каждый столбец
        int _t = order[new_bas_idx];
        order[new_bas_idx] = order[m_prev];
        order[m_prev] = _t;


        // std::cout << "tut5" << std::endl;
        canone(A_can, B_can);

        // std::cout << "tut6" << std::endl;
        ss.A = A_can; ss.B = B_can; ss.order = order; ss.c_ext = c_ext;
    }
}

SystemState push_variable(
    SystemState ss) {

    int m_curr = ss.A.rows();
    int n_1_curr = ss.A.cols();
    int n_curr = n_1_curr - m_curr;

    auto x_curr = ss.x_curr;

    // finding indmax{x_i}
    int max_fract_idx = -1;
    for (int i = 0; i < m_curr; i++) {
        //std::cout << (x_curr[i] - floor(x_curr[i])) << ";;";
        if (ss.order[i] >= ss.opt_vars) {
            continue;
        }
        if (max_fract_idx < 0 || (x_curr[ss.order[i]] - floor(x_curr[ss.order[i]]))
        > (x_curr[ss.order[max_fract_idx]] - floor(x_curr[ss.order[max_fract_idx]]))) {
            max_fract_idx = i;
        }
    }

    //std::cout << max_fract_idx << "=>" << ss.order[max_fract_idx] << std::endl;
    //std::cout << x_curr[ss.order[max_fract_idx]] - floor(x_curr[ss.order[max_fract_idx]]) << std::endl;
    Matrix A = Matrix(m_curr + 1, n_1_curr + 1);
    Matrix B = Matrix(m_curr + 1, 1);

    for (int i = 0; i < m_curr; i++) {
        for (int j = 0; j < n_1_curr; j++) {
            A[i][j] = ss.A[i][j];
        }
        A[i][n_1_curr] = 0;
        B[i][0] = ss.B[i][0];
    }
    B[m_curr][0] = x_curr[ss.order[max_fract_idx]] - floor(x_curr[ss.order[max_fract_idx]]);
    
    for (int j = 0; j < m_curr; j++) {
        A[m_curr][j] = 0;
    }

    // std::cout << "tut612" << std::endl;
    for (int j = m_curr; j < n_1_curr; j++) {
        A[m_curr][j] = (ss.A[max_fract_idx][j] - floor(ss.A[max_fract_idx][j]));
    }
    A[m_curr][n_1_curr] = -1;

    // std::cout << "tut613" << std::endl;
    std::vector<int> order(ss.order.size() + 1, ss.order.size());
    for (int i = 0; i < ss.order.size(); i++) {
        order[i] = ss.order[i];
    }

    // std::cout << "tut614" << std::endl;
    std::vector<long double> c_ext(ss.c_ext.size() + 1, 0);
    for (int i = 0; i < ss.c_ext.size(); i++) {
        c_ext[i] = ss.c_ext[i];
    }

    return SystemState(NOT_YET, order, A, B, c_ext, ss.opt_vars);
}

SystemState simplex_init(
    const Matrix& A, const Matrix& B,
    const std::vector<long double>& c) {
    // step 0
    int n = A.cols();
    int m = A.rows();
    int n_1 = n + m;

    // Создаём доп переменные, чтобы перейти от общей к канон форме
    Matrix E = Matrix(m, m);
    for (int i = 0; i < m; i++) {
        E[i][i] = 1;
    }
    auto A_can = concat(E, A);

    // step 1
    // Расчёт каноничной формы (E|A'|B')
    auto B_can = B;
    canone(A_can, B_can);
    std::vector<long double> c_ext(n_1, 0.);
    for (int i = 0; i < c.size(); i++) {
        c_ext[m + i] = c[i];
    }

    // step 2
    // запоминаем какой переменной соответсвует каждый столбец
    std::vector<int> order(n_1);
    for (int i = 0; i < n_1; i++) {
        order[(n + i) % n_1] = i;
    }
    return SystemState(NOT_YET, order, A_can, B_can, c_ext, c.size());
}

SystemState simplex_step(SystemState ss) {
    auto A_can = ss.A;
    auto B_can = ss.B;
    auto order = ss.order;
    auto c_ext = ss.c_ext;
    
    int n_1 = A_can.cols();
    int m = A_can.rows();
    int n = n_1 - m;

    std::vector<long double> x(n_1, 0.0);
    for (int i = 0; i < B_can.rows(); i++) {
        x[order[i]] = B_can[i][0];
    }
    std::vector<long double> z;

    int result = NOT_YET;
    int max_iterations = 20;

    for (int iter = 0; iter < max_iterations; iter++) {
        std::cout << "##" << iter << std::endl;
        // step 3
        // считаем данные в симплекс-таблице
        z.clear();
        z.resize(n_1, 0.);
        std::vector<long double> Delta(n_1, 0.);
        std::vector<long double> t(m, 0.);
        calc_table(A_can, c_ext, z, Delta);

        // step 4
        // проверяем получено ли решение
        int solve_row_idx, solve_col_idx;
        result = analyze(A_can, B_can, x, c_ext, z, Delta, t, solve_row_idx, solve_col_idx);
        //print_table(A_can, B_can, c_ext, z, Delta, t, order);

        if (result != NOT_YET) {
            // сохраняем решение к двойственной задаче.
            // они лежат в z_i, соответсвующих небазисным переменным x_i
            /*for (int i = 0; i < order.size(); i++) {
                if (order[i] >= n)
                    x[order[i]] = z[i];
            }*/
            return SystemState(result, order, A_can, B_can, c_ext, ss.opt_vars, x);
        }
        //std::cout << "A:" << A_can;
        //std::cout << "B:" << B_can;

        // step 5
        // переходим к новому базису, обновляем таблицу, приводим к канон виду
        swap_cols(A_can, c_ext, solve_row_idx, solve_col_idx);
        canone(A_can, B_can);
        // запоминаем какой переменной соответсвует каждый столбец
        int _t = order[solve_col_idx];
        order[solve_col_idx] = order[solve_row_idx];
        order[solve_row_idx] = _t;


        // пересчитываем x (с учетом смены столбцов)
        x.clear();
        x.resize(n_1, 0.0);
        // небазисные зануляются, остальные равны правой части B
        for (int i = 0; i < B_can.rows(); i++) {
            x[order[i]] = B_can[i][0];
        }
    }
    //std::cout << "ITERATION LIMIT!" << std::endl;
    return SystemState(ITER_LIMIT, order, A_can, B_can, c_ext, ss.opt_vars, x);
}

// Не понадобилось, т.к. я по своему храню порядок базисных и небаз переменных.
// Я заменил эту функцию -> swap_cols() + canone()
void update_ab(Matrix& A, Matrix& B, int solve_row_idx, int solve_col_idx) {
    auto A_1 = A;
    auto B_1 = B;

    for (int i = 0; i < A.rows(); i++) {
        if (i == solve_row_idx)
            continue;
        B_1[i][0] = B[i][0] - B[solve_row_idx][0] * A[i][solve_col_idx] / A[solve_row_idx][solve_col_idx];
        for (int j = 0; j < A.cols(); j++) {
            A_1[i][j] = A[i][j] - A[solve_row_idx][j] * A[i][solve_col_idx] / A[solve_row_idx][solve_col_idx];
        }
    }

    A_1[solve_row_idx][solve_col_idx] = 1;
    for (int i = 0; i < A.rows(); i++) {
        if (i == solve_row_idx)
            continue;
        A_1[i][solve_col_idx] = 0;
    }

    B_1[solve_row_idx][0] = B[solve_row_idx][0] / A[solve_row_idx][solve_col_idx];
    for (int j = 0; j < A.cols(); j++) {
        A_1[solve_row_idx][j] = A[solve_row_idx][j] / A[solve_row_idx][solve_col_idx];
    }
    A = A_1;
    B = B_1;
    //std::cout << "A_1:" << A_1;
    //std::cout << "B_1:" << B_1;
}

void swap_cols(Matrix& A, std::vector<long double>& c, int p, int q)  {
    long double _t = c[p];
    c[p] = c[q];
    c[q] = _t;

    for (int i = 0; i < A.rows(); i++) {
        _t = A[i][p];
        A[i][p] = A[i][q];
        A[i][q] = _t;
    }
}

int analyze(const Matrix& A, const Matrix& B, const std::vector<long double>& x,
    const std::vector<long double>& c, const std::vector<long double>& z,
    const std::vector<long double>& Delta, std::vector<long double>& t,
    int& solve_row_idx, int& solve_col_idx) {
    bool all_pos = true;
    solve_col_idx = 0;
    for (int i = 0; i < Delta.size(); i++) {
        if (Delta[i] >= 0) {
            continue;
        }
        all_pos = false;
        if (x[i] < 0) {
            return INFINITE_SOL;
        }
        if (Delta[i] < Delta[solve_col_idx]) {
            solve_col_idx = i;
        }
    }
    if (all_pos) {
        return FINITE_SOL;
    }
    solve_row_idx = -1;
    for (int i = 0; i < t.size(); i++) {
        if (A[i][solve_col_idx] <= 0) {
            t[i] = -1;
            continue;
        }
        t[i] = B[i][0] / A[i][solve_col_idx];
        if (solve_row_idx < 0 || (t[i] < t[solve_row_idx] && t[i] > 0)) {
            solve_row_idx = i;
        }
    }
    //std::cout << solve_row_idx << "," << solve_col_idx << ":" << A[solve_row_idx][solve_col_idx] << std::endl;

    return NOT_YET;
}

void print_table(const Matrix& A, const Matrix& B, const std::vector<long double>& c,
    const std::vector<long double>& z, const std::vector<long double>& Delta,
    const std::vector<long double>& t, const std::vector<int>& order) {
    // использую для красоты: << std::setw(10) << std::setprecision(4)
    std::cout << std::setfill('-') << std::setw(c.size()*9+32) << "" << std::endl << std::setfill(' ');
    std::cout << std::setw(5) << "" << std::setw(9) << std::left << "c" << std::setw(9) << "";
    for (int i = 0; i < c.size(); i++) {
        std::cout << std::setw(9) << std::setprecision(4) << c[i];
    }
    std::cout << std::setw(9) << std::left << "t" << std::endl;
    std::cout << std::endl;
    std::cout << std::setw(5) << "basis" << std::setw(9) << "" << std::setw(9) << std::left << "b";
    for (int i = 0; i < A.cols(); i++) {
        std::cout << "a^" << std::setw(7) << std::left << order[i] + 1;
    }
    std::cout << std::endl;
    for (int i = 0; i < A.rows(); i++) {
        std::cout << "a_" << std::setw(3) << std::left << order[i] + 1
            << std::setw(9) << std::setprecision(4) << c[i]
            << std::setw(9) << std::setprecision(4) << B[i][0];
        for (int j = 0; j < A.cols(); j++) {
            std::cout << std::setw(9) << std::setprecision(4) << A[i][j];
        }
        std::cout << std::setw(9) << std::setprecision(4) << std::left
            //<< (t[i] > 0 ? std::to_string(t[i]) : "") << std::endl;
            << t[i] << std::endl;
    }
    std::cout << std::setw(5) << "z" << std::setw(9) << "" << std::setw(9) << "?";
    for (int i = 0; i < z.size(); i++) {
        std::cout << std::setw(9) << std::setprecision(4) << z[i];
    }
    std::cout << std::endl;
    std::cout << std::setw(5) << "Delta" << std::setw(9) << "" << std::setw(9) << "";
    for (int i = 0; i < Delta.size(); i++) {
        std::cout << std::setw(9) << std::setprecision(4) << Delta[i];
    }
    std::cout << std::endl;
    std::cout << std::setfill('-') << std::setw(c.size() * 9 + 32) << "" << std::endl << std::setfill(' ');
}

void calc_table(const Matrix& A, const std::vector<long double>& c,
    std::vector<long double>& z, std::vector<long double>& Delta) {
    for (int i = 0; i < A.cols(); i++) {
        z[i] = 0;
        for (int j = 0; j < A.rows(); j++) {
            z[i] += c[j] * A[j][i];
        }
        Delta[i] = z[i] - c[i];
    }
}

Matrix concat(const Matrix& A, const Matrix& B) {
    Matrix C = Matrix(A.rows(), A.cols() + B.cols());
    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < A.cols(); j++) {
            C[i][j] = A[i][j];
        }
        for (int j = A.cols(); j < A.cols() + B.cols(); j++) {
            C[i][j] = B[i][j - A.cols()];
        }
    }
    return C;
}

Matrix canone(Matrix& A, Matrix& B)
{
    //assert A.rows() == B.rows();
    int n = A.cols();
    int m = A.rows();
    //assert n >= m;

    Matrix A_new = concat(A, B);
    //std::cout << "Прямой" << std::endl;
    //Прямой ход (Зануление нижнего левого угла)
    for (int k = 0; k < m; k++) //k-номер строки
    {
        // std::cout << "tut0" << std::endl;
        int not_zero_row = k;
        //while (not_zero_row < m and A_new[not_zero_row][k] == 0.) {
        while (A_new[not_zero_row][k] == 0.) {
            not_zero_row++;
        }
        //assert not_zero_row < m
        //swap k_row and not_zero_row
        // std::cout << "tut1" << std::endl;
        if (k < not_zero_row) {
            auto _t = A_new[k];
            A_new[k] = A_new[not_zero_row];
            A_new[not_zero_row] = _t;


            _t = A[k];
            A[k] = A[not_zero_row];
            A[not_zero_row] = _t;
        }

        // std::cout << "tut2" << std::endl;
        for (int i = 0; i < n + 1; i++) //i-номер столбца
            A_new[k][i] = A_new[k][i] / A[k][k]; //Деление k-строки на первый член !=0 для преобразования его в единицу

        // std::cout << "tut3" << std::endl;
        for (int i = k + 1; i < m; i++) //i-номер следующей строки после k
        {
            double K = A_new[i][k] / A_new[k][k]; //Коэффициент
            for (int j = 0; j < n + 1; j++) //j-номер столбца следующей строки после k
                A_new[i][j] = A_new[i][j] - A_new[k][j] * K; //Зануление элементов матрицы ниже первого члена, преобразованного в единицу
        }

        // std::cout << "tut4" << std::endl;
        for (int i = 0; i < m; i++) //Обновление, внесение изменений в начальную матрицу
        {
            for (int j = 0; j < n; j++) {
                A[i][j] = A_new[i][j];
            }
            B[i][0] = A_new[i][n];
        }
    }

    //std::cout << "Обр" << std::endl;
    //Обратный ход (Зануление верхнего правого угла)
    for (int k = m - 1; k > -1; k--) //k-номер строки
    {
        for (int i = n; i > -1; i--) //i-номер столбца
            A_new[k][i] = A_new[k][i] / A[k][k];
        for (int i = k - 1; i > -1; i--) //i-номер следующей строки после k
        {
            double K = A_new[i][k] / A_new[k][k];
            for (int j = n; j > -1; j--) //j-номер столбца следующей строки после k
                A_new[i][j] = A_new[i][j] - A_new[k][j] * K;
        }
    }

    for (int i = 0; i < m; i++) //Обновление, внесение изменений в начальную матрицу
    {
        for (int j = 0; j < n; j++) {
            A[i][j] = A_new[i][j];
        }
        B[i][0] = A_new[i][n];
    }

    return A_new;
};