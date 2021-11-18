#include <Rcpp.h>
using namespace Rcpp;

//naming conventions: m - mean, s - standard deviation, n - count (= event_length * 3012)

// [[Rcpp::export]]
double sd_c(double x_m, double x_s, double x_n,double y_m, double y_s, int y_n)
{
	double al, var, tmp_sd;
            al=x_n+y_n;

            tmp_sd=al*((x_n-1)*(x_s*x_s)+(y_n-1)*(y_s*y_s))+y_n*x_n*(x_m-y_m)*(x_m-y_m);
            var=tmp_sd/(al*(al-1));
            
            return(sqrt(var));
}

// [[Rcpp::export]]
double sd_d(NumericVector sd,NumericVector mean, NumericVector freq)
{
	double i , x1, nx, ny, sx, sy, mx, my;
	x1 = mean.size();

	if(x1==1){
		sx=sd[0];
		return(sx);
	}
	mx=mean[0];
	sx=sd[0];
	nx=freq[0];
	for (i =1;i<x1;i++)
	{
		my=mean[i];
		sy=sd[i];
		ny=freq[i];

		sx=sd_c(mx,sx,nx,my,sy,ny);
		mx=(mx*nx+my*ny)/(nx+ny);
		nx=nx+ny;
	}

	return(sx);
}	

// [[Rcpp::export]]
//This function aggregates the event level mean for a k-mer segment, if multiple are available, on the same read 
double mean_c(NumericVector mean, NumericVector freq)
{
	double i , x1, nx, ny, mx, my;
	//obtain size of input vector 'mean'
	x1 = mean.size(); 

	if(x1==1){
		nx=mean[0];
		return(nx);
	}
	mx=mean[0];
	nx=freq[0];
	for (i =1;i<x1;i++)
	{
		my=mean[i];
		ny=freq[i];

		mx=(mx*nx+my*ny)/(nx+ny);
		nx=nx+ny;
	}

	return(mx);
}	

//This code converts ASCII characters into integers
// [[Rcpp::export]]
void myascii(NumericVector qual, NumericVector qualn) {
	for (int i = 0; i<qual.size(); i++) {qualn[i] = (int)qual[i]-1;}
}

//This code generates match/mismatch and deletion columns
// [[Rcpp::export]]
void mymatmis(CharacterVector op, CharacterVector ref, CharacterVector base, NumericVector mat, NumericVector mis, NumericVector del) {
	for (int i=0;i < mat.size();i++) {
		if (op[i] == "M") {
			if (ref[i] == base[i]) {mat[i] = 1;mis[i] = 0;del[i] = 0;}
			else {mat[i] = 0;mis[i] = 1;del[i] = 0;}
		}
		else if (op[i] == "D") {mat[i] = 0;mis[i] = 0;del[i] = 1;}
		else {mat[i] = 0;mis[i] = 0;del[i] = 0;}
	}
}


//This code generates insertion column postive strand
// [[Rcpp::export]]
void insertfunP(CharacterVector op, NumericVector inse) {
  for (int i=0;i < inse.size();i++) {
    if (op[i] == "I") {
      int k=1;
      while (op[i+k]=="I"){
        k++;
      }
      inse[i-1]=k;
    }}
}

// [[Rcpp::export]]
void countnb5(CharacterVector refe, CharacterVector strand, NumericVector pos, NumericVector countld1, NumericVector countlg1,
              NumericVector countld2, NumericVector countlg2) {
  for (int i=0;i < pos.size();i++) {
    if ((refe[i] == "A" && strand[i] == "p")|(refe[i] == "T" && strand[i] == "n")) {
      int k=0;
      while (pos[i+k+1]==pos[i]+1){
        k++;
      }
      countld1[i]=k;
      
      int h=0;
      while (pos[i-h-1]==pos[i]-1){
        h++;
      }
      countlg1[i]=h;
      
      int m=0;
      while (pos[i+k+m+1]==pos[i]+2){
        m++;
      }
      countld2[i]=m;
      
      int n=0;
      while (pos[i-h-n-1]==pos[i]-2){
        n++;
      }
      countlg2[i]=n;
    }
  }
}
    
//This code generates insertion column negative strand
// [[Rcpp::export]]
void insertfunN(CharacterVector op, NumericVector inse) {
  for (int i=0;i < inse.size();i++) {
    if (op[i] == "I") {
      int k=1;
      while (op[i-k]=="I"){
        k++;
      }
      inse[i+1]=k;
    }}
}


// [[Rcpp::export]]
  std::string rev (std::string dna) {
    int n = dna.size();
    std::string result(n, 'x');
    for (int i = 0; i < n; i++) {
      result[n - 1 - i]  = dna[i];
    }
    return result;
}

// [[Rcpp::export]]
std::string complement (std::string dna) {
  int n = dna.size();
  std::string result(n, 'x');
  for (int i = 0; i < n; i++) {
    switch (dna[i]) {
    case 'A' :
      result[i] = 'T';
      break;
    case 'T' :
      result[i] = 'A';
      break;
    case 'G' :
      result[i] = 'C';
      break;
    case 'C' :
      result[i] = 'G';
      break;
    }
  }
  return result;
}

// [[Rcpp::export]]
 void complement_char(CharacterVector base, CharacterVector base1) {
  for (int i=0;i < base.size();i++) {
      if (base[i] == "A") {base1[i] = "T";}
      else if (base[i] == "T") {base1[i] = "A";}
      else if (base[i] == "C") {base1[i] = "G";}
      else if (base[i] == "G") {base1[i] = "C";}
      else {base1[i]="na";}
 }
 }
 
 // [[Rcpp::export]]
Rcpp::StringVector StringV_to_vstrings(Rcpp::StringVector StringV){
    std::vector<std::string> vstrings(StringV.size());
    int i;
    for (i = 0; i < StringV.size(); i++){
        vstrings[i] = StringV(i);
    }
    
    Rcpp::StringVector StringV2(StringV.size());
    StringV2 = vstrings;
    return(StringV2);
}


// [[Rcpp::export]]
void rev_comp(StringVector kmer, StringVector kmer1){
    std::vector<std::string> kmerstr(kmer.size());
    for (int i = 0; i < kmer.size(); i++){
        kmerstr[i] = kmer(i);
    }
    
    for (int i=0; i< kmer.size(); i++){
        kmer1[i]=rev(complement(kmerstr[i]));
    }
}

// [[Rcpp::export]]
NumericVector indexgroup_big(NumericVector ind, int sizevector){
  int bin=10;
  float b=(sizevector/10.0);
  int nr = floor(b);
  for(int r=0; r<(bin);r++){
    for(int i=(r*nr); i<((r+1)*nr); i++){
      ind[i]=r+1;
    }
  }
  for(int k=((bin-1)*nr); k<sizevector; k++){
    ind[k]=bin;
  }
  return ind;
}

// [[Rcpp::export]]
NumericVector indexgroup10(NumericVector ind, int sizevector){
  int bin;
  float b=(sizevector/10.0);
  bin = floor(b);
  for(int r=0; r<(bin);r++){
    for(int i=(r*10); i<((r+1)*10); i++){
        ind[i]=r+1;
    }
  }
  for(int k=((bin-1)*10); k<sizevector; k++){
      ind[k]=bin;
  }
  return ind;
}



// [[Rcpp::export]]
NumericMatrix transform_each_pos5 (NumericVector x, NumericMatrix y, int ncol) {
  int nx = x.size();
  int ny = y.nrow();
  // copy x to 2nd row
  //for (int i = 0; i < nx; i++) result[nx + i] = x[i];
  // 
  NumericVector resultb1;
  NumericVector resulta1;
  NumericVector resultb2;
  NumericVector resulta2;
  for (int i = 0; i < ny; i++) {
    if ((y(i,0)==x[0]) && (y(i,1)==(x[1]-1))) {
      resultb1.push_back(i);}
    if ((y(i,0)==x[0]) && (y(i,1)==(x[1]+1))) {
      resulta1.push_back(i);}
    if ((y(i,0)==x[0]) && (y(i,1)==(x[1]-2))) {
      resultb2.push_back(i);}
    if ((y(i,0)==x[0]) && (y(i,1)==(x[1]+2))) {
      resulta2.push_back(i);}
  }
  
  int nb1=resultb1.size();
  int na1=resulta1.size();
  int nb2=resultb2.size();
  int na2=resulta2.size();
  NumericMatrix result(5*nb1*nb2*na1*na2,ncol);
  for (int i = 0; i < nb1*nb2*na1*na2; i++) {
    result(5*i+2,_) = x;
  }

  for (int i = 0; i < nb2; i++){
    for (int j = 0; j < nb1; j++){
      for (int h = 0; h < na1; h++){
        for (int k = 0; k < na2; k++){
          int ta = (nb1*na1*na2)*i+(na1*na2)*j+na2*h+k;
          result(5*ta,_) = y(resultb2[i],_);
          result(5*ta+1,_) = y(resultb1[j],_);
          result(5*ta+3,_) = y(resulta1[h],_);
          result(5*ta+4,_) = y(resulta2[k],_);
          }
      }
    }
  }
  
  
  return result;
}

// [[Rcpp::export]]
NumericVector mysample(NumericVector v,int n) { 
  return sample(v, n, true); 
}


// [[Rcpp::export]]
NumericMatrix transform_each_pos5_regular (NumericVector x, NumericMatrix y, int ncol) {
  int ny = y.nrow();
  // copy x to 2nd row
  //for (int i = 0; i < nx; i++) result[nx + i] = x[i];
  //
  NumericVector resultb1;
  NumericVector resulta1;
  NumericVector resultm;
  NumericVector resultb2;
  NumericVector resulta2;

  for (int i = 0; i < ny; i++) {
    if ((y(i,0)==x[0]) && (y(i,1)==(x[1]-1))) {
      resultb1.push_back(i);}
    if ((y(i,0)==x[0]) && (y(i,1)==(x[1]+1))) {
      resulta1.push_back(i);}
    if ((y(i,0)==x[0]) && (y(i,1)==(x[1]))) {
      resultm.push_back(i);}
    if ((y(i,0)==x[0]) && (y(i,1)==(x[1]-2))) {
      resultb2.push_back(i);}
    if ((y(i,0)==x[0]) && (y(i,1)==(x[1]+2))) {
      resulta2.push_back(i);}
  }
  
  NumericVector b2;
  NumericVector b1;
  NumericVector m;
  NumericVector a1;
  NumericVector a2;
  
  NumericMatrix result(50,ncol);
  for (int i = 0; i < 10; i++) {
    b2=mysample(resultb2,1);
    b1=mysample(resultb1,1);
    m=mysample(resultm,1);
    a1=mysample(resulta1,1);
    a2=mysample(resulta2,1);
    result(5*i,_)=y(b2[0],_);
    result(5*i+1,_)=y(b1[0],_);
    result(5*i+2,_) = y(m[0],_);
    result(5*i+3,_)=y(a1[0],_);
    result(5*i+4,_)=y(a2[0],_);
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix transform_each_pos5_big (NumericVector x, NumericMatrix y,int ncol) {
  int ny = y.nrow();
  // copy x to 2nd row
  //for (int i = 0; i < nx; i++) result[nx + i] = x[i];
  //
  NumericVector resultb1;
  NumericVector resulta1;
  NumericVector resultm;
  NumericVector resultb2;
  NumericVector resulta2;
  
  for (int i = 0; i < ny; i++) {
    if ((y(i,0)==x[0]) && (y(i,1)==(x[1]-1))) {
      resultb1.push_back(i);}
    if ((y(i,0)==x[0]) && (y(i,1)==(x[1]+1))) {
      resulta1.push_back(i);}
    if ((y(i,0)==x[0]) && (y(i,1)==(x[1]))) {
      resultm.push_back(i);}
    if ((y(i,0)==x[0]) && (y(i,1)==(x[1]-2))) {
      resultb2.push_back(i);}
    if ((y(i,0)==x[0]) && (y(i,1)==(x[1]+2))) {
      resulta2.push_back(i);}
  }
  
  NumericVector b2;
  NumericVector b1;
  NumericVector a1;
  NumericVector a2;
  
  NumericMatrix result(50,ncol);
  for (int i = 0; i < 10; i++) {
    b2=mysample(resultb2,1);
    b1=mysample(resultb1,1);
    a1=mysample(resulta1,1);
    a2=mysample(resulta2,1);
    result(5*i,_)=y(b2[0],_);
    result(5*i+1,_)=y(b1[0],_);
    result(5*i+2,_) = y(resultm[i],_);
    result(5*i+3,_)=y(a1[0],_);
    result(5*i+4,_)=y(a2[0],_);
  }
  return result;
}

// [[Rcpp::export]]
CharacterVector replace_list(StringVector pat, CharacterVector repl, StringVector x){
  int nx=x.size();
  CharacterVector y(nx);
  for (int i=0; i<nx; i++){
    for (int j=0; j<pat.size(); j++){
      if (x[i]==pat[j]){y[i]=repl[j];}
    }
  }
  return y;
}



