/*********************************************
 * Concert Model
 * Autor: Charles Paulino
 * Data de criação: 04-06-2019 
 * Problem - MODELO FLEET-ICt ATUAL
 *********************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <ilcplex/ilocplex.h>
using namespace std;

typedef IloArray<IloNumVarArray> NumVarMatrix2;
typedef IloArray<NumVarMatrix2> NumVarMatrix3;
typedef IloArray<NumVarMatrix3> NumVarMatrix4;
typedef vector<vector<int> > matrix2;
typedef vector<matrix2> matrix3;

class Mono{
    public:
        /// Model parameter
        int SS; // tamanho da matriz dist_bases
        int n; //nós de demanda
        int m; //bases potenciais
        int v; //períodos
        int Q; // número máximo de bases USA's + USB's
        int w; //tipos de ambulâncias (2)
        vector<int> P; // qt de cada tipo de ambulância
        vector<int> S; // tempo de resposta de cada tipo de ambulância
        vector<int> C; // capacidade de cada base (2 veículos)
        matrix2 d; // distância entre j e i/
        matrix3 q; // demanda de cada ponto i, por cada tipo de ambulância u, em cada período t

        /// Model parameter 
        IloEnv env; // create a common environment (build + solve model)
        IloModel mod; // create LP, QP, MILP, MIQP model (with variables, constraints, etc.)
        IloCplex cplex; //Solve model with CPLEX
        IloObjective fo; //função objetivo
        IloNumVarArray z;
        NumVarMatrix3 y;
        NumVarMatrix3 x;

    public:
        Mono(){
            mod = IloModel(env);
            cplex = IloCplex(mod);
        }
        void Read_data(char name[]);
        void Create_model();
        void Display_configuration();
};

void  Mono::Read_data(char name[]){
    ifstream arq(name);
    if (!arq.is_open()){
        cout << "Error openning file: " << name << endl;
        arq.close();
        exit(EXIT_FAILURE);
    }
	arq >> SS; //
    arq >> n;
    arq >> m;
    arq >> v;
    arq >> Q;
    arq >> w;
    P = vector<int>(w);
    S = vector<int>(w);
    C = vector<int>(m);
    
    d = matrix2(m);
    for(int j = 0; j < m; j++){
		d[j] = vector<int>(n);
	}
	
    q = matrix3(n);
    for(int i = 0; i < n; i++){
		q[i] = matrix2(w);
		for(int u = 0; u < w; u++){
			q[i][u] = vector<int>(v);
        }}   
        

    for (int u = 0; u < w; u++){
		arq >> P[u];
	}
    for (int u = 0; u < w; u++){
		arq >> S[u];
    }
    for (int j = 0; j < m; j++){
		arq >> C[j];
    }
    
	int h;
	int l;
	for (int i = 0; i <= SS; i++){ 
			arq >> h
				>> l;
			arq >> d[h][l];
		}

    for (int i = 0; i < n; i++){
		for (int t = 0; t < v; t++){
			for (int u = 0; u < w; u++){ 
			arq >> q[i][u][t];
		}}}    
  	arq.close();
}

void Mono::Create_model(){
    z = IloNumVarArray(env, m, 0.0, 1.0, ILOINT);
    y = NumVarMatrix3(env, n);
    for(int i = 0; i < n; i++){
		y[i] = NumVarMatrix2(env, w);
		for(int u = 0; u < w; u++){
			y[i][u] = IloNumVarArray(env, v, 0.0, 1.0, ILOINT);
	}}
    x = NumVarMatrix3(env, m);
    for(int j = 0; j < m; j++){
		x[j] = NumVarMatrix2(env, w);
		for(int u = 0; u < w; u++){
			x[j][u] = IloNumVarArray(env, v, 0.0, 1.0, ILOINT);
	}} 	

    // ==============================
    // objective function
    // ==============================
    
    IloExpr expfo(env);
    for (int i = 0; i <  n; i++){
		for (int u = 0; u < w; u++){
			for (int t = 0; t < v; t++){
        expfo += q[i][u][t] * y[i][u][t]; 
    }}}
    IloAdd(mod, IloMaximize(env, expfo));
    
    //===================================
    // constraint: forall(u in U, t in T) sum(j in J) x[j][u][t] <= P[u]
    //===================================
        for (int u = 0; u < w; u++){
			for (int t = 0; t < v; t++){
			IloExpr r1(env);
			for (int j = 0; j < m; j++){
			r1 += x[j][u][t];
			}
            mod.add(r1 <= P[u]);
            r1.end();
        }}
            
     //====================================================
    // constraint forall(j in J, t in T) sum(u in U) x[j][u][t] <= C[j] * z[j]
     //====================================================
        for (int j = 0; j < m; j++){
			for (int t = 0; t < v; t++){	   
        IloExpr r2(env);
        for (int u = 0; u < w; u++){
			r2 += x[j][u][t];
        }
        mod.add(r2 <= C[j] * z[j]);
        r2.end();
    }}

     //====================================================
    // constraint sum(j in J) z[j] <= Q
     //====================================================   
        IloExpr r3(env);
        for (int j = 0; j < m; j++){
			r3 += z[j];
        }
        mod.add(r3 <= Q);
        r3.end();
    
     //====================================================
    // constraint forall(j in J, t in T, u in U) x[j][u][t] <= z[j]; 
     //====================================================
        for (int j = 0; j < m; j++){
			for (int t = 0; t < v; t++){
				for (int u = 0; u < w; u++){
					mod.add(x[j][u][t] <= z[j]);
    }}}

     //====================================================
    // constraint forall(i in I, u in U, t in T) y[i][u][t] <= sum (j in J) x[j][u][t];
     //====================================================
		for (int i = 0; i < n; i++){	
			for (int u = 0; u < w; u++){
				for (int t = 0; t < v; t++){
					IloExpr r5(env);
					for (int j = 0; j < m; j++){
						if(d[j][i] <= S[u]){
							r5 += x[j][u][t];
							}}
							mod.add(r5 >= y[i][u][t]);
							r5.end();
    }}}

}
void Mono::Display_configuration(){
    printf("BASES:");
    for(int j = 0; j < m; j++){
        if( cplex.getValue(z[j]) >= 0.5){
            cout<<j + 1<<"\t";
        }
    }
    cout<<endl;
   printf("\n ALOCAÇÃO DAS AMBULÂNCIAS: \n");
    for(int j = 0; j < m; j++){
		for(int u = 0; u < w; u++){
			for(int t = 0; t < v; t++){
				 if( cplex.getValue(x[j][u][t]) >= 0.5) {
				cout<<j + 1<<"\t"<<u + 1<<"\t"<<t + 1<<endl;
            }}}}
    
    	printf("\n TOTAL DE BASES INSTALADAS: \n");
			int soma3 = 0;
			for(int j = 0; j < m; j++){
				if( cplex.getValue(z[j]) >= 0.5){
					soma3 += cplex.getValue(z[j]);
					}}
					cout<<soma3<<endl;
 
		printf("\n DEMANDA COBERTA PARA CADA TIPO DE AMBULÂNCIA EM CADA NÍVEL: \n");
				for(int u = 0; u < w; u++){
					int soma4 = 0;
					int soma5 = 0;
					for(int i = 0; i < n; i++){
						for(int t = 0; t < v; t++){
						soma4 += q[i][u][t];
						if( cplex.getValue(y[i][u][t]) >= 0.5){
						soma5 += q[i][u][t];
						}}}
						cout<<soma4<<"\t"<<soma5<<endl;
					}
}

int main (int argc, char *argv[]){

    Mono *mono=new Mono();

    try{

        // Read_main_arg(argc, argv);
        cout<<"Openning file "<<argv[1]<<endl;
        if(argc > 1){
			mono->Read_data(argv[1]);
		}
		else{
			mono->Read_data("Data_FLEET-ICt_2t");
		}

        mono->Create_model();

        /// ==========================
        /// configurações do cplex
        /// ==========================
	//	mono->cplex.setParam(IloCplex::EpGap, 0.5); // Definindo uma tolerancia de GAP
        mono->cplex.setWarning(mono->env.getNullStream());
   //     mono->cplex.setOut(mono->env.getNullStream()); //se comentar imprime os logs

        ///==============================
        /// Resolvendo o problema
        ///==============================

        IloTimer crono(mono->env);
        crono.start();
        mono->cplex.solve();
        crono.stop();

		cout<<"Status "<<mono->cplex.getStatus()<<endl; /// Imprimindo o status
		mono->Display_configuration(); /// Imprimindo a configuração
    
        ///=====================================
        /// Salvando os resultados no arquivo
        ///=====================================

        FILE *fp;
        fp = fopen("Results.txt","aw+");
        fprintf(fp,"%s\t%d\t%d\t%1.2f\t%1.6f\n",  argv[1], mono->n, mono->m,   mono->cplex.getObjValue(), (double) crono.getTime());
        printf("%s\t%d\t%d\t%1.2f\t%1.6f\n",  argv[1], mono->n, mono->m,   mono->cplex.getObjValue(), (double) crono.getTime());
    } 
    catch (IloException& ex) {
        cerr << "Error: " << ex << endl;
    }
    return 0;
}
