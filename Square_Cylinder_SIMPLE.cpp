#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
using namespace std;
void TDMA(int y,int i,int m,int n,int cx,int cy, double**a,double**b,double**c,double**d,double**p_star)
{
    int j;
    double P[y+2],Q[y+2];

    for(j=1; j<y+1; j++)
    {
        if( i>m-cx-1 && j>n-cy-1 && j<n+cy+1 && i<m+cx+1)
        {
            continue;
        }
        if (j==1)
        {
            P[j]=b[i][j]/a[i][j];
            Q[j]=d[i][j]/a[i][j];
        }
        else
            P[j]=b[i][j]/(a[i][j]-(c[i][j]*P[j-1]));
        Q[j]=(d[i][j]+(c[i][j]*Q[j-1]))/(a[i][j]-(c[i][j]*P[j-1]));

    }
    p_star[i][y]=Q[y];
    for(j=y-1; j>=1; j--)
    {
        if( i>m-cx-1 && j>n-cy-1 && j<n+cy+1 && i<m+cx+1)
        {
            continue;
        }
        p_star[i][j]=Q[j]+(P[j]*p_star[i+1][j]);
    }
}
//Function to check whether value are getting calculated outside the cylinder
int square(int i, int j, int s,int m,int cx,int cy, int n)
{
    if (s == 0)   // u-velocity control volume
    {
        if (i>=(m-cx) && i<(m+cx+2) && j>=(n-cy) && j<(n+cy+1))
        {
            return 1;
        }
        else
            return 0;
    }
    if (s == 1)   // v-velocity control volume
    {
        if (i>=(m-cx) && i<(m+cx+1) && j>=(n-cy) && j<(n+cy+2))
        {
            return 1;
        }
        else
            return 0;
    }
    if (s == 2)   // pressure  control volume
    {
        if (i>=(m-cx) && i<(m+cx+1) && j>=(n-cy) && j<(n+cy+1))
        {
            return 1;
        }
        else
            return 0;
    }
}
int main()
{
    int i,j,itr,x,y,n,m,cx,cy,k;
    double Re,mu,rho,length,height,dy,dx,c_side,U,alpha_u,alpha_v,alpha_p,p_sum,res=1,sum;
    c_side=0.01; //side length of cylinder
    height=10*c_side;
    length=60*c_side;
    y=61;//number of cells in y direction
    x=y*6;//number of cells in x direction
    n=y/2+1;//cell position of square center
    m=x/4;// square center in x direction
    dx=length/x;
    dy=height/y;
    cx=c_side/(2*dx);//no of cells from center to cylinder side in x
    cy=c_side/(2*dy);// no of cells from center to cylinder side in y
cout<<"SIMPLE algorithm for staggered grid"<<endl;
cout<<"Flow over square cylinder in the duct"<<endl;
cout<<"Grid size 366*61"<<endl;

    U=1;
    mu=0.1;
    cout<<"Inlet velocity = "<<U<<"m/s dynamic viscosity = "<<mu<<endl;
cout<<"Select density (1,5,10,20,50) for Reynold number (0.1,0.5,1,2,5) respectively"<<endl;
cin>>rho;
    Re=(U*c_side*rho)/mu;
    alpha_u=0.9; // under-relaxation factor for u momentum equation
    alpha_v=0.9; // under-relaxation factor for v momentum equation
    alpha_p=0.001; // under-relaxation factor for pressure equation

    double p[x+2][y+2],u[x+2][y+2],v[x+2][y+2],u_star[x+2][y+2],v_star[x+2][y+2];
    double apu[x+2][y+2],apv[x+2][y+2],app[x+2][y+2],bp[x+2][y+2],ae[x+2][y+2],aw[x+2][y+2],an[x+2][y+2],as[x+2][y+2];
    double fe[x+2][y+2],fw[x+2][y+2],fn[x+2][y+2],fs[x+2][y+2],df[x+2][y+2],pre[x+2][y+2],de,dw,ds,dn;
    //diffusion coefficients
    de=(mu*dy)/dx;
    dw=(mu*dy)/dx;
    dn=(mu*dx)/dy;
    ds=(mu*dx)/dy;
    cout<<"Reynold Number =  "<<Re<<endl;
    // Dynamic memory allocation
    double** p_star = new double*[x+2];
    double** a = new double*[x+2];
    double** b = new double*[x+2];
    double** c = new double*[x+2];
    double** d = new double*[x+2];

    for (i=0; i<x+2 ; i++)
    {
        p_star[i] = new double[y+2];
        a[i] = new double[y+2];
        b[i] = new double[y+2];
        c[i] = new double[y+2];
        d[i] = new double[y+2];
    }

    for(i=0; i<(x+2); i++)
    {
        for(j=0; j<(y+2); j++)
        {
            u[i][j]= 1;
            u_star[i][j]=1; // initial guess for u velocity
            v[i][j]= 0;
            v_star[i][j]=0; // initial guess for v velocity
            p[i][j]= 1;     // initial guess for pressure
            p_star[i][j]=0;
            pre[i][j]=0;
            apu[i][j]= 0;
            apv[i][j]= 0;
            app[i][j]= 0;
            an[i][j]= 0;
            as[i][j]= 0;
            ae[i][j]= 0;
            aw[i][j]= 0;
            bp[i][j]= 0;
            fn[i][j]= 0;
            fs[i][j]= 0;
            fe[i][j]= 0;
            fw[i][j]= 0;
            df[i][j]= 0;
        }
    }
// Uniform velocity at inlet
    for(j=1; j<y+1; j++)
    {
        u[1][j]=U;
        u[x+1][j]=u[x][j];
        u_star[1][j]=U;
        u_star[x+1][j]=u_star[x][j];
        if(j>1)
        {
            v[x+1][j]=v[x][j];
            v_star[x+1][j]=v[x][j];
        }
    }
    for(i=m-cx; i<(m+cx+2); i++)
    {
        for(j=n-cy; j<(n+cy+1); j++)
        {
            u[i][j]=0;
            u_star[i][j]=0;          // square cylinder interior
        }
    }
    for(i=m-cx; i<(m+cx+1); i++)
    {
        for(j=m-cy; j<(m+cy+2); j++)
        {
            v[i][j]=0;
            v_star[i][j]=0;     // square cylinder interior
        }
    }

    itr=0;
    cout<<setprecision(6)<<fixed<<endl;
    //*********SIMPLE ALOGIRITHM************//
    while(abs(res)>0.0001)
    {
        for (i=2; i<x+1; i++)
        {
            for (j=1; j<y+1; j++)
            {
                //Mass fluxes calculation
                fn[i][j] = 0.5*dx*rho*(v_star[i-1][j+1] + v_star[i][j+1]);
                fs[i][j] = 0.5*dx*rho*(v_star[i-1][j] + v_star[i][j]);
                fe[i][j] = 0.5*dy*rho*(u_star[i][j] + u_star[i+1][j]);
                fw[i][j] = 0.5*dy*rho*(u_star[i-1][j] + u_star[i][j]);
                df[i][j] = fn[i][j] - fs[i][j] + fe[i][j] - fw[i][j];
                //Hybrid Scheme implementation
                an[i][j] = max(-fn[i][j],max((dn-fn[i][j]/2),0.00));
                as[i][j] = max(fs[i][j],max(ds+fs[i][j]/2,0.00));
                ae[i][j] = max(-fe[i][j],max(de-fe[i][j]/2,0.00));
                aw[i][j] = max(fw[i][j],max(dw+fw[i][j]/2,0.00));
            }
        }
        for(i=2; i<x+1; i++)
        {
            for(j=1; j<y+1; j++)
            {
                bp[i][j]=dy*(p[i-1][j]-p[i][j]); //source term
                apu[i][j]= an[i][j] + as[i][j] + ae[i][j] + aw[i][j] + df[i][j];
            }
        }
        // Gauss Seidel method to solve u momentum equation
        for (k=0; k<3000; k++)
        {
            for (i=2; i<(x+1); i++)
            {
                for (j=1; j<(y+1); j++)
                {
                    if (square(i,j,0,m,cx,cy,n)==0)
                    {
                        u[i][j] = (1-alpha_u)*u_star[i][j] + (alpha_u/apu[i][j])*(an[i][j]*u[i][j+1] + as[i][j]*u[i][j-1] + ae[i][j]*u[i+1][j] + aw[i][j]*u[i-1][j] + bp[i][j]);
                        if(isnan(u[i][j]))
                        {
                            cout<<i<<j<<"t";
                            return 0;
                        }
                    }
                }
            }
            //Boundary conditions
            for(i=2; i<x+1; i++)
            {
                u[i][y+1]=-u[i][y]; //duct top wall
                u[i][0]=-u[i][1]; //duct bottom wall
            }
            for(j=1; j<y+1; j++)
            {
                u[1][j]=U;      // inlet
                u[x+1][j]=u[x][j]; // zero gradient outlet
            }
            for(i=m-cx; i<=(m+cx+1); i++)
            {
                u[i][n-cy]=-u[i][n-cy-1]; // square cylinder bottom wall
                u[i][n+cy]=-u[i][n+cy+1];    // square cylinder top wall
            }
        }
        // solve y-momentum equation

        // discretised coefficients
        for (i=1; i<x+1; i++)
        {
            for (j=2; j<y+1; j++)
            {
                //Mass flux calculation
                fn[i][j] = 0.5*dx*rho*(v_star[i][j] + v_star[i][j+1]);
                fs[i][j] = 0.5*dx*rho*(v_star[i][j-1] + v_star[i][j]);
                fe[i][j] = 0.5*dy*rho*(u_star[i+1][j-1] + u_star[i+1][j]);
                fw[i][j] = 0.5*dy*rho*(u_star[i][j-1] + u_star[i][j]);
                df[i][j] = fn[i][j] - fs[i][j] + fe[i][j] - fw[i][j];
                //Hybrid Scheme implementation
                an[i][j] = max(-fn[i][j], max(dn - fn[i][j]/2, 0.00));
                as[i][j] = max(fs[i][j], max(ds + fs[i][j]/2, 0.00));
                ae[i][j] = max(-fe[i][j], max(de - fe[i][j]/2, 0.00));
                aw[i][j] = max(fw[i][j], max(dw + fw[i][j]/2, 0.00));
            }
        }
        for(i=1; i<x+1; i++)
        {
            for(j=2; j<y+1; j++)
            {
                bp[i][j]=dx*(p[i][j-1]-p[i][j]);//source term
                apv[i][j]= an[i][j] + as[i][j] + ae[i][j] + aw[i][j] + df[i][j];
            }
        }
         // Gauss Seidel method to solve u momentum equation
        for (k=0; k<3000; k++)
        {
            for (i=1; i<x+1; i++)
            {
                for (j=2; j<y+1; j++)
                {
                    if  (square(i,j,1,m,cx,cy,n)==0)
                    {
                        v[i][j] = (1-alpha_v)*v_star[i][j] + (alpha_v/apv[i][j])*(an[i][j]*v[i][j+1] + as[i][j]*v[i][j-1] + ae[i][j]*v[i+1][j] + aw[i][j]*v[i-1][j] + bp[i][j]);
                        if(isnan(v[i][j]))
                        {
                            cout<<i<<" "<<j<<"\t";
                            return 0;
                        }
                    }
                }
            }
             //Boundary conditions
            for(i=1; i<x+1; i++)
            {
                v[i][y+1]=0; //duct top wall
                v[i][1]=0;  //duct bottom wall
            }
            for(j=2; j<y+1; j++)
            {
                v[x+1][j]=v[x][j]; //outlet
                v[0][j]=-v[1][j]; //inlet
            }
            for(j=n-cy; j<=n+cy+1; j++)
            {
                v[m-cx][j]=-v[m-cx-1][j]; //cylinder front wall
                v[m+cx][j]=-v[m+cx+1][j]; //cylinder back wall
            }
        }
        // solve pressure correction equation

        // discretised coefficients
        for(i=1; i<x; i++)
        {
            for(j=1; j<y; j++)
            {
                an[i][j]=(dx*dx*alpha_v)/apv[i][j+1];
                as[i][j+1]=(dx*dx*alpha_v)/apv[i][j+1];
            }
            an[i][y]=0; //top wall
            as[i][1]=0; // bottom wall
        }
        for(i=1; i<x; i++)
        {
            for(j=1; j<y+1; j++)
            {
                ae[i][j]=(dy*dy*alpha_u)/apu[i+1][j];//check for alpha
                aw[i+1][j]=(dy*dy*alpha_u)/apu[i+1][j];
                aw[1][j]=0; //inlet
            }
        }
        for(i=m-cx; i<m+cx+1; i++)
        {
            an[i][n-cy-1]=0; //bottom wall of cylinder
            as[i][n+cy+1]=0; //top wall of cylinder
        }
        for(j=n-cy; j<n+cy+1; j++)
        {
            ae[m-cx-1][j]=0; //front wall of cylinder
            aw[m+cx+1][j]=0; //back wall of cylinder
        }

        // iterative solver
        for(i=1; i<x+1; i++)
        {
            for(j=1; j<y+1; j++)
            {
                bp[i][j]=(v[i][j]-v[i][j+1])*dx+(u[i][j]-u[i+1][j])*dy;
                app[i][j]=an[i][j]+as[i][j]+ae[i][j]+aw[i][j];
                p_star[i][j]=0;
            }
        }
       //Line by Line TDMA solver for pressure correction equation
        for (k=0; k<3000 ; k++)
        {
            for (i=1; i<x; i++)
            {
                for (j=1; j<y+1; j++)
                {
                    p_star[x][j]=0; //fixed outlet pressure
                    if (square(i,j,2,m,cx,cy,n)==0)
                    {
                        if(i==m-cx-1 && j>=n-cy && j<=n+cy)
                        {
                            d[i][j]=bp[i][j]+(p_star[i-1][j]*aw[i][j]);
                            a[i][j]=an[i][j]+as[i][j]+aw[i][j];
                            b[i][j]=an[i][j];
                            c[i][j]=as[i][j];
                        }
                        if(i==m+cx+1 && j>=n-cy && j<=n+cy)
                        {
                            d[i][j]=bp[i][j]+(p_star[i+1][j]*ae[i][j]);
                            a[i][j]=an[i][j]+as[i][j]+ae[i][j];
                            b[i][j]=an[i][j];
                            c[i][j]=as[i][j];
                        }
                        if(j==n-cy-1 && i>=m-cx && i<=m+cx)
                        {
                            d[i][j]=bp[i][j]+(p_star[i-1][j]*aw[i][j])+(p_star[i+1][j]*ae[i][j]);
                            a[i][j]=as[i][j]+ae[i][j]+aw[i][j];
                            b[i][j]=0;
                            c[i][j]=as[i][j];
                        }
                        if(j==n+cy+1 && i>=m-cx && i<=m+cx)
                        {
                            d[i][j]=bp[i][j]+(p_star[i-1][j]*aw[i][j])+(p_star[i+1][j]*ae[i][j]);
                            a[i][j]=an[i][j]+ae[i][j]+aw[i][j];
                            b[i][j]=an[i][j];
                            c[i][j]=0;
                        }
                        if(j==1 && i==1)
                        {
                            d[i][j]=bp[i][j]+(p_star[i+1][j]*ae[i][j]);
                            a[i][j]=ae[i][j]+an[i][j];
                            b[i][j]=an[i][j];
                            c[i][j]=0;
                        }
                        if(j==1 && i==x-1)
                        {
                            d[i][j]=bp[i][j]+(p_star[i-1][j]*aw[i][j]);
                            a[i][j]=aw[i][j]+an[i][j]+ae[i][j];
                            b[i][j]=an[i][j];
                            c[i][j]=0;
                        }
                        if(i==1 && j==y)
                        {
                            d[i][j]=bp[i][j]+(p_star[i+1][j]*ae[i][j]);
                            a[i][j]=ae[i][j]+as[i][j];
                            b[i][j]=0;
                            c[i][j]=as[i][j];
                        }
                        if(i==x-1 && j==y)
                        {
                            d[i][j]=bp[i][j]+(p_star[i-1][j]*aw[i][j]);
                            a[i][j]=aw[i][j]+as[i][j]+ae[i][j];
                            b[i][j]=0;
                            c[i][j]=as[i][j];
                        }
                        if(i==1 && 1<j && j<y)
                        {
                            d[i][j]=bp[i][j]+(p_star[i+1][j]*ae[i][j]);
                            a[i][j]=an[i][j]+as[i][j]+ae[i][j];
                            b[i][j]=an[i][j];
                            c[i][j]=as[i][j];
                        }
                        if(i==x-1 && 1<j && j<y)
                        {
                            d[i][j]=bp[i][j]+(p_star[i-1][j]*aw[i][j]);
                            a[i][j]=an[i][j]+as[i][j]+aw[i][j];
                            b[i][j]=an[i][j];
                            c[i][j]=as[i][j];
                        }
                        if(j==1 && 1<i && i<x-1)
                        {
                            d[i][j]=bp[i][j]+(p_star[i-1][j]*aw[i][j])+(p_star[i+1][j]*ae[i][j]);
                            a[i][j]=an[i][j]+ae[i][j]+aw[i][j];
                            b[i][j]=an[i][j];
                            c[i][j]=0;
                        }
                        if(1<i && i<x-1 && j==y)
                        {
                            d[i][j]=bp[i][j]+(p_star[i-1][j]*aw[i][j])+(p_star[i+1][j]*ae[i][j]);
                            a[i][j]=as[i][j]+ae[i][j]+aw[i][j];
                            b[i][j]=0;
                            c[i][j]=as[i][j];
                        }
                        else
                            d[i][j]=bp[i][j]+(p_star[i-1][j]*aw[i][j])+(p_star[i+1][j]*ae[i][j]);
                        a[i][j]=as[i][j]+ae[i][j]+aw[i][j]+an[i][j];
                        b[i][j]=an[i][j];
                        c[i][j]=as[i][j];
                    }
                }
                TDMA(y,i,m,n,cx,cy,a,b,c,d,p_star); //function to perform TDMA solver
            }
        }

        // pressure and velocity correction
        for(i=1; i<x+1; i++)
        {
            for(j=1; j<y+1; j++)
            {
                if (square(i,j,2,m,cx,cy,n)==0)
                {
                    p[i][j]=p[i][j]+alpha_p*p_star[i][j];
                }
                if(j>1)
                {
                    if (square(i,j,1,m,cx,cy,n)==0)
                    {
                        v[i][j]=v[i][j]+(dx/apv[i][j])*(p_star[i][j-1]-p_star[i][j]);
                        v_star[i][j]=v[i][j];
                    }
                }
            }
            if(i>1)
            {
                for(j=1; j<y+1; j++)
                {
                    if (square(i,j,0,m,cx,cy,n)==0)
                    {
                        u[i][j]=u[i][j]+(dy/apu[i][j])*(p_star[i-1][j]-p_star[i][j]);
                        u_star[i][j]=u[i][j];
                    }
                }
            }
        }

        // mass residual
        sum=0,res=0;
        for(i=1; i<x+1; i++)
        {
            for(j=1; j<y+1; j++)
            {
                sum=(v[i][j+1]-v[i][j])*dx+(u[i+1][j]-u[i][j])*dy;
                res=res+sum;
            }
        }

        cout<<itr<<"\t"<<"mass residual"<<'\t'<<res<<endl;
        itr=itr+1;
    }
 //Post processing
    for(i=m-cx; i<m+cx+2; i++)
    {
        for(j=n-cy; j<n+cy+1; j++)
        {
            u[i][j]=0;

        }
        if(i<m+cx+2)
        {
            for(j=n-cy; j<n+cy+2; j++)
            {
                v[i][j]=0;
            }
        }
    }
// Drag coefficient calculation
    p_sum=0;
    for(j=n-cy-1; j<=n+cy+1; j++)
    {
        p_sum=p_sum+(p[m-cx-1][j]-p[m+cx+1][j]);
    }
    cout<<"Drag Coefficient = "<<(p_sum*dy)/(0.5*rho*U*U*c_side)<<endl;
   // Result file for plotting
    ofstream output;
    output.open("result.csv");
    output<<"x,y,z,p,u,v"<<endl;

    for(j=1; j<y+1; j++)
    {
        for(i=1; i<x+1; i++)
        {
            output<<i-1<<","<<j-1<<","<<"0,"<<p[i][j]<<","<<u[i][j]<<","<<v[i][j]<<endl;

        }
    }
    return 0;
}


