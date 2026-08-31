#include <string>
#include <stdio.h>
#include <amgx_c.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;



// Function to generate different configurations based on the provided name
const char *getConfig(const std::string &configName)
{
    char *cfg;

    // Generate JSON string based on the configName using if-else statements
    if (configName == "AMG_CLASSICAL_PMIS")
    {
        cfg = "{\"config_version\":2,\"determinism_flag\":1,\"solver\":{\"print_grid_stats\":1,\"store_res_history\":1,\"obtain_timings\":1,\"solver\":\"GMRES\",\"print_solve_stats\":1,\"preconditioner\":{\"interpolator\":\"D2\",\"solver\":\"AMG\",\"cycle\":\"V\",\"smoother\":{\"relaxation_factor\":1,\"scope\":\"jacobi\",\"solver\":\"JACOBI_L1\"},\"presweeps\":2,\"postsweeps\":2,\"selector\":\"PMIS\",\"coarsest_sweeps\":2,\"coarse_solver\":\"NOSOLVER\",\"max_iters\":1,\"max_row_sum\":0.9,\"min_coarse_rows\":2,\"scope\":\"amg_solver\",\"max_levels\":24,\"print_grid_stats\":1,\"aggressive_levels\":1,\"interp_max_elements\":4},\"max_iters\":100,\"monitor_residual\":1,\"gmres_n_restart\":10,\"convergence\":\"RELATIVE_INI_CORE\",\"tolerance\":1e-06,\"norm\":\"L2\"}}";

    }
    else if (configName == "GMRES_AMG_D2")
    {
        cfg = "{\"config_version\":2,\"determinism_flag\":1,\"exception_handling\":1,\"solver\":{\"print_grid_stats\":1,\"store_res_history\":1,\"solver\":\"GMRES\",\"print_solve_stats\":1,\"obtain_timings\":1,\"preconditioner\":{\"interpolator\":\"D2\",\"print_grid_stats\":1,\"solver\":\"AMG\",\"smoother\":\"JACOBI_L1\",\"presweeps\":2,\"selector\":\"PMIS\",\"coarsest_sweeps\":2,\"coarse_solver\":\"NOSOLVER\",\"max_iters\":1,\"interp_max_elements\":4,\"min_coarse_rows\":2,\"scope\":\"amg_solver\",\"max_levels\":24,\"cycle\":\"V\",\"postsweeps\":2},\"max_iters\":200,\"monitor_residual\":1,\"gmres_n_restart\":10,\"convergence\":\"ABSOLUTE\",\"tolerance\":1e-05,\"norm\":\"L2\"}}";

    }
    else if (configName == "PBICGSTAB_NOPREC")
    {
        cfg = "{\"config_version\":2,\"solver\":{\"preconditioner\":{\"scope\":\"amg_solver\",\"solver\":\"NOSOLVER\"},\"use_scalar_norm\":1,\"solver\":\"PBICGSTAB\",\"print_solve_stats\":1,\"obtain_timings\":1,\"monitor_residual\":1,\"convergence\":\"ABSOLUTE\",\"scope\":\"main\",\"tolerance\":1e-5,\"norm\":\"L2\"}}";

    }
    else if (configName == "V-cheby-smoother")
    {
        cfg = "{\"config_version\":2,\"determinism_flag\":1,\"solver\":{\"scope\":\"main\",\"print_grid_stats\":1,\"solver\":\"AMG\",\"scaling\":\"DIAGONAL_SYMMETRIC\",\"interpolator\":\"D2\",\"aggressive_levels\":0,\"interp_max_elements\":4,\"coarse_solver\":\"NOSOLVER\",\"print_solve_stats\":1,\"obtain_timings\":1,\"max_iters\":400,\"monitor_residual\":1,\"convergence\":\"ABSOLUTE\",\"max_levels\":100,\"cycle\":\"V\",\"smoother\":{\"solver\":\"CHEBYSHEV\",\"preconditioner\":{\"solver\":\"NOSOLVER\",\"max_iters\":1},\"max_iters\":1,\"chebyshev_polynomial_order\":4,\"chebyshev_lambda_estimate_mode\":2},\"tolerance\":1e-06,\"norm\":\"L2\",\"presweeps\":0,\"postsweeps\":1}}";

    }
    else
    {
        std::cerr << "Error: Unknown configuration name." << std::endl;
        return nullptr;
    }

    return cfg;
}

void readCSRMatrix(AMGX_matrix_handle matrix, AMGX_vector_handle rhs, AMGX_vector_handle soln, std::string problem)
{
    std::string path = "problems/" + problem + "/";
    std::string filename = path + "csr_matrix.txt";
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }

    int nnz = 0;
    int n = 0;
    std::string line;

    // Read ia: (number of values)
    std::getline(file, line);
    line = line.substr(3);
    n = std::stoi(line);
    n -= 1;

    // Allocate memory for the arrays
    int *row_ptrs = new int[n + 1];

    // Read ia: (integer values separated by spaces)
    std::getline(file, line);
    std::istringstream iss(line);
    for (int i = 0; i <= n; i++)
    {
        iss >> row_ptrs[i];
        row_ptrs[i] -= 1;
    }

    // Read ja: (integer values separated by spaces)
    std::getline(file, line);
    line = line.substr(3);
    nnz = std::stoi(line);
    cout << nnz << endl;
    int *col_indices = new int[nnz];
    double *values = new double[nnz];

    iss.clear();
    std::getline(file, line);
    iss.str(line);

    for (int i = 0; i < nnz; i++)
    {
        iss >> col_indices[i];
        col_indices[i] -= 1;
    }

    // Read values: (double values separated by spaces)
    iss.clear();
    std::getline(file, line);
    std::getline(file, line);
    iss.str(line);

    for (int i = 0; i < nnz; i++)
    {
        iss >> values[i];
    }

    file.close();

    AMGX_matrix_upload_all(matrix, n, nnz, 1, 1, row_ptrs, col_indices, values, NULL);

    delete[] row_ptrs;
    delete[] col_indices;
    delete[] values;

    // Read rhs: (double values separated by newlines)
    double *rhs_values = new double[n];
    double *soln_values = new double[n];

    filename = path + "rhs.txt";
    std::ifstream rhs_file(filename);
    if (!rhs_file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }

    for (int i = 0; i < n; i++)
    {
        rhs_file >> rhs_values[i];
    }

    rhs_file.close();

    // Soln is initialized to 0
    for (int i = 0; i < n; i++)
    {
        soln_values[i] = 0.0;
    }

    AMGX_vector_upload(rhs, n, 1, rhs_values);
    AMGX_vector_upload(soln, n, 1, soln_values);
}


void readSolution(double *solution, int size, std::string problem)
{
    string filename = "problems/" + problem + "/solution.txt";
    exit(0);
    std::ifstream solution_file("solution.txt");

    if (!solution_file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }

    for (int i = 0; i < size - 1; i++)
    {
        solution_file >> solution[i];
    }

    solution_file.close();
}

void gerarMatrizS(AMGX_matrix_handle matrix, AMGX_vector_handle rhs, AMGX_vector_handle soln, double stepSize)
{
    double a = -1, b = 1; // intervalo

    double alfa = 2.0;
    double beta = 2.0;

    int numeroPontos = ((b - a) / stepSize) + 1;
    cout << "numeroPontos: " << numeroPontos << endl;

    int n = numeroPontos - 2;

    double h = stepSize * stepSize;

    int nnz = 3 * n - 2; // elementos não-nulos - n+(n-1)+(n-1) = 3n-2

    int *row_ptrs = new int[n + 1];
    int *col_indices = new int[nnz];
    double *values = new double[nnz];
    double *rhs_values = new double[n];
    double *soln_values = new double[n]();

    row_ptrs[0] = 0;
    int index = 0;

    for (int i = 0; i < n; i++)
    {
        if (i > 0)
        {
            col_indices[index] = i - 1;
            values[index] = -1.0;
            index++;
        }

        col_indices[index] = i;
        values[index] = 2.0;
        index++;

        if (i < n - 1)
        {
            col_indices[index] = i + 1;
            values[index] = -1.0;
            index++;
        }

        row_ptrs[i + 1] = index;
    }

    double vd = -h * 2.0 + alfa;

    rhs_values[0] = vd;
    for (int i = 1; i < n - 1; i++)
    {
        rhs_values[i] = -2.0 * h;
    }
    rhs_values[n - 1] = vd;

    AMGX_matrix_upload_all(matrix, n, nnz, 1, 1, row_ptrs, col_indices, values, NULL);
    AMGX_vector_upload(rhs, n, 1, rhs_values);
    AMGX_vector_upload(soln, n, 1, soln_values);

    delete[] row_ptrs;
    delete[] col_indices;
    delete[] values;
    delete[] rhs_values;
    delete[] soln_values;
}

void gerarMatrizArquivo(double stepSize)
{

    double a = -1, b = 1; // intervalo

    double alfa = 2.0;
    double beta = 2.0;

    int numeroPontos = ((b - a) / stepSize) + 1;
    cout << "numeroPontos: " << numeroPontos << endl;

    /*
        numeroPontos = 9 são 7 valores + 2 valores de contorno
        na resolução do sistema não uso o contorno então faço -2
    */

    int n = numeroPontos - 2;

    // double stepSize = (b - a) / (numeroPontos - 1); // valor de cada incremento, 9 pontos mas 8 intervalos
    // cout << "n: " << n << endl;
    // cout << "STEP SIZE: " << stepSize << endl;

    double h = stepSize * stepSize;
    // cout << "h: " << h << endl;

    int nnz = 3 * n - 2; // elementos não-nulos - n+(n-1)+(n-1) = 3n-2
    // cout << "NNZ: " << nnz << endl;
    int countI = 1, countF = 3;

    ofstream arquivoMtx;

    arquivoMtx.open("../examples/matrix3.mtx");
    arquivoMtx << fixed << scientific << setprecision(15);
    arquivoMtx << "%%MatrixMarket matrix coordinate real general symmetric" << endl;
    arquivoMtx << "%%AMGX 1 1 sorted rhs" << endl;
    arquivoMtx << n << " " << n << " " << nnz << endl;
    arquivoMtx << 1 << " " << 1 << " " << 2.0 << endl;
    arquivoMtx << 1 << " " << 2 << " " << -1.0 << endl;

    for (int i = 2; i <= n - 1; i++)
    {
        for (int j = countI; j <= countF; j++)
        {
            if (i == j)
            {
                arquivoMtx << i << " " << j << " " << 2.0 << endl;
            }
            else
            {
                arquivoMtx << i << " " << j << " " << -1.0 << endl;
            }
        }
        countI++;
        countF++;
    }

    arquivoMtx << n << " " << n - 1 << " " << -1.0 << endl;
    arquivoMtx << n << " " << n << " " << 2.0 << endl;

    double vd = -h * 2.0 + alfa;
    // cout << "RHS: " << vd << endl;

    arquivoMtx << vd << endl;
    for (int i = 0; i < n - 2; i++)
    {
        arquivoMtx << -2.0 * h << endl;
    }
    arquivoMtx << vd;

    arquivoMtx.close();
}

void calcular(const char **argv, double stepSize)
{
    /*
        0 -> nome do programa
        1 -> -c
        2 -> arquivo config
        3 -> -s
        4 -> stepSize
    */

    // cout << "Config: " << argv[2] << endl;

    // cout << "STEP SIZE: " << stepSize << endl;

    AMGX_initialize();

    // gerarMatrizArquivo(stepSize);

    // All of the objects are initialized using a default Resources.
    AMGX_matrix_handle matrix;
    AMGX_vector_handle rhs;
    AMGX_vector_handle soln;
    AMGX_resources_handle rsrc;
    AMGX_solver_handle solver;
    AMGX_config_handle config;

    // arquivo de configuração passado como parametro: ../src/configs/FGMRES_AGGREGATION_JACOBI.json
    AMGX_config_create_from_file(&config, argv[2]);
    AMGX_resources_create_simple(&rsrc, config);

    AMGX_matrix_create(&matrix, rsrc, AMGX_mode_dDDI);
    AMGX_vector_create(&rhs, rsrc, AMGX_mode_dDDI);
    AMGX_vector_create(&soln, rsrc, AMGX_mode_dDDI);

    double *data = new double[100000000];

    /*
        1 d -> the first letter h or d specifies whether the matrix data (and subsequent linear solver algorithms) will run on the host or device
        2 D(or F) -> specifies the precision (double or float) of the Matrix data
        3 D The third D or F specifies the precision (double or float) of any Vector (including right-hand side or unknown vectors)
        4 I The last I specifies that 32-bit int types are used for all indices
    */
    AMGX_solver_create(&solver, rsrc, AMGX_mode_dDDI, config);

    // Next, data is uploaded from the application (or set, in the case of the solution vector which is initialized to all zeroes).
    //  If these are not specified than rhs=[1,...,1]^T and (initial guess) sol=[0,...,0]^T.
    gerarMatrizS(matrix, rhs, soln, stepSize);
    // AMGX_read_system(matrix, rhs, soln, "../examples/matrix3.mtx");
    //readCSRMatrix(matrix, rhs, soln);

    // AMGX_write_system(matrix, rhs, soln, "./output.system.mtx");

    AMGX_solver_setup(solver, matrix);

    AMGX_solver_solve_with_0_initial_guess(solver, rhs, soln);

    AMGX_vector_download(soln, data);

    int sol_size, sol_bsize;
    AMGX_vector_get_size(soln, &sol_size, &sol_bsize);

    double *solution = new double[sol_size];
    // readSolution(solution, sol_size);

    ofstream plotSol;

    plotSol.open("plotSol.csv");
    // plotSol << 2 << endl;
    double errAcc = 0;
    for (int i = 0; i < sol_size; ++i)
    {
        // printf("%f \n",data[i]);
        // plotSol << data[i]<<","<<solution[i]<<","<<pow(solution[i]- data[i],2)<< endl;
        plotSol << data[i] << endl;
        errAcc += pow(solution[i] - data[i], 2);
    }
    // plotSol << "Erro: "<<errAcc << endl;
    //  plotSol << 2 << endl;
    plotSol.close();

    AMGX_solver_destroy(solver);
    AMGX_vector_destroy(soln);
    AMGX_vector_destroy(rhs);
    AMGX_matrix_destroy(matrix);
    AMGX_resources_destroy(rsrc);

    AMGX_SAFE_CALL(AMGX_config_destroy(config));
    AMGX_SAFE_CALL(AMGX_finalize());

    delete[] data;
}

void calcular(const char **argv, string problem)
{

    /*
        0 -> nome do programa
        1 -> -c
        2 -> config
        3 -> -p
        4 -> problem
    */

    // cout << "Config: " << argv[2] << endl;

    // cout << "Problem: " << problem << endl;

    AMGX_initialize();

    // All of the objects are initialized using a default Resources.
    AMGX_matrix_handle matrix;
    AMGX_vector_handle rhs;
    AMGX_vector_handle soln;
    AMGX_resources_handle rsrc;
    AMGX_solver_handle solver;
    AMGX_config_handle config;

    // arquivo de configuração passado como parametro: ../src/configs/FGMRES_AGGREGATION_JACOBI.json
    std::string configString = "PBICGSTAB_NOPREC";
    const char *cfg = getConfig(configString);
    AMGX_config_create(&config, cfg);
    cout << "Created config " << configString << endl;
    AMGX_resources_create_simple(&rsrc, config);

    AMGX_matrix_create(&matrix, rsrc, AMGX_mode_dDDI);
    AMGX_vector_create(&rhs, rsrc, AMGX_mode_dDDI);
    AMGX_vector_create(&soln, rsrc, AMGX_mode_dDDI);

    double *data = new double[100000000];

    /*
        1 d -> the first letter h or d specifies whether the matrix data (and subsequent linear solver algorithms) will run on the host or device
        2 D(or F) -> specifies the precision (double or float) of the Matrix data
        3 D The third D or F specifies the precision (double or float) of any Vector (including right-hand side or unknown vectors)
        4 I The last I specifies that 32-bit int types are used for all indices
    */
    AMGX_solver_create(&solver, rsrc, AMGX_mode_dDDI, config);

    // Next, data is uploaded from the application (or set, in the case of the solution vector which is initialized to all zeroes).
    //  If these are not specified than rhs=[1,...,1]^T and (initial guess) sol=[0,...,0]^T.
    // gerarMatrizS(matrix, rhs, soln, stepSize);
    // AMGX_read_system(matrix, rhs, soln, "../examples/matrix3.mtx");
    readCSRMatrix(matrix, rhs, soln, problem);

    // AMGX_write_system(matrix, rhs, soln, "./output.system.mtx");

    AMGX_solver_setup(solver, matrix);

    AMGX_solver_solve_with_0_initial_guess(solver, rhs, soln);

    AMGX_vector_download(soln, data);

    int sol_size, sol_bsize;
    AMGX_vector_get_size(soln, &sol_size, &sol_bsize);

    double *solution = new double[sol_size];
    // readSolution(solution, sol_size, problem);

    ofstream plotSol;

    plotSol.open("plotSol.csv");
    // plotSol << 2 << endl;
    double errAcc = 0;
    for (int i = 0; i < sol_size; ++i)
    {
        // printf("%f \n",data[i]);
        // plotSol << data[i]<<","<<solution[i]<<","<<pow(solution[i]- data[i],2)<< endl;
        plotSol << data[i] << endl;
        errAcc += pow(solution[i] - data[i], 2);
    }
    // plotSol << "Erro: "<<errAcc << endl;
    //  plotSol << 2 << endl;
    plotSol.close();

    AMGX_solver_destroy(solver);
    AMGX_vector_destroy(soln);
    AMGX_vector_destroy(rhs);
    AMGX_matrix_destroy(matrix);
    AMGX_resources_destroy(rsrc);

    AMGX_SAFE_CALL(AMGX_config_destroy(config));
    AMGX_SAFE_CALL(AMGX_finalize());

    delete[] data;
}

int main(int argc, const char **argv)
{
    // exemplo de chamada: examples/meu_teste2 -c GMRES_AMG_D2 -p cube-poisson

    /*
        0 -> nome do programa
        1 -> -p
        2 -> problem
    */

    // cout << "Config: " << argv[2] << endl;
    int numeroPontos;
    double stepSize;
    istringstream ss(argv[2]);
    if (std::string(argv[1]) == "-n")
    {
        ss >> numeroPontos;
        // cout << "Número Pontos: " << numeroPontos << endl;
        // gerarMatrizN(numeroPontos);
    }
    else if (std::string(argv[1]) == "-s")
    {
        ss >> stepSize;
        cout << "Step Size: " << stepSize << endl;
        // gerarMatrizS(stepSize);
        calcular(argv, stepSize);
    }
    else if (std::string(argv[1]) == "-p")
    {
        std::string problem = argv[2];
        cout << "Solving problem: " << problem << endl;
        calcular(argv, problem);
    }

    return 0;
}
