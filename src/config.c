#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "color_mdff.h"
#include "config.h"
#include "constants.h"
#include "cell.h"
#include "tools.h"
#include "io.h"
#include "global.h"
#include "pbc.h"
#include "verlet.h"
#include "kinetic.h"
#include "field.h"

void pre_alloc_config(){
    for (int it=0;it<NTYPEMAX;it++)
    {
        atypit[it]= malloc((MAX_LEN+1)*sizeof (char*));
        massit[it]= 1.0;
        nionit[it]= 0;
    }
}

//-----------------------------------------------------------------------------
void alloc_config(){

    onenion=1.0/( (double) nion) ;
    onethirdnion=1.0/( 3.0*(double) nion) ;
    massia= malloc(nion*sizeof *massia );
    invemassia= malloc(nion*sizeof *invemassia );
    atypia= malloc(nion*sizeof *atypia);
    if (atypia==NULL){
        io_node printf("**atype malloc failure\n");
    }
    for (int ia=0;ia<nion;ia++)  {
        atypia[ia]=NULL;
        atypia[ia]=malloc((MAX_LEN+1)*sizeof(*atypia));
        if (atypia[ia] == NULL){
            io_node printf("atype %d malloc failure\n",ia);
        }
    }
    vx= malloc(nion*sizeof(*vx));
    vy= malloc(nion*sizeof(*vy));
    vz= malloc(nion*sizeof(*vz));
    rx= malloc(nion*sizeof(*rx));
    ry= malloc(nion*sizeof(*ry));
    rz= malloc(nion*sizeof(*rz));
    fx= malloc(nion*sizeof(*fx));
    fy= malloc(nion*sizeof(*fy));
    fz= malloc(nion*sizeof(*fz));
    rxs= malloc(nion*sizeof(*rxs));
    rys= malloc(nion*sizeof(*rys));
    rzs= malloc(nion*sizeof(*rzs));
    for(int ia=0; ia<nion; ia++) {
        massia[ia]=1.0;
        invemassia[ia]=1.0;
        rx[ia]=0.0;ry[ia]=0.0;rz[ia]=0.0;
        vx[ia]=0.0;vy[ia]=0.0;vz[ia]=0.0;
        fx[ia]=0.0;fy[ia]=0.0;fz[ia]=0.0;
        rxs[ia]=0.0;rys[ia]=0.0;rzs[ia]=0.0;
    }
    if (lverletL) { 
        verlet_nb=allocate_verletlist("vnlnb");
        verlet_nb->cut=cutshortrange+skindiff; 
    }

    if (lcoul) {
        qia=malloc(nion*sizeof(*qia));
        dipia=malloc(nion*sizeof(*dipia));
        quadia=malloc(nion*sizeof(*quadia));
        if (lverletL) {
            verlet_coul=allocate_verletlist("vnlcoul");
            verlet_coul->cut=cutlongrange+skindiff; 
        }
    }

}

//-----------------------------------------------------------------------------
void free_config(){
    for (int it=0;it<ntype;it++)  {
        free(atypit[it]);
    }
    free(massia);
    free(invemassia);
    for (int ia=0;ia<nion;ia++)  {
        free(atypia[ia]);
    }
    free(atypia);
    free(rx);free(ry);free(rz);
    free(vx);free(vy);free(vz);
    free(fx);free(fy);free(fz);
    free(rxs);free(rys);free(rzs);
    if (lverletL) free_verletlist("vnlnb");

    if (lcoul) {
        free(qia);
        free(dipia);
        free(quadia);
    }

}

//-----------------------------------------------------------------------------
bool config_input_ok(){
    
    int sumnionit=sum(nionit,ntype);
    if ( sumnionit != nion) {
        pError("nionit and nion does not correspond\n");
        printf("%d %d\n",sumnionit,nion);
        return false;
    }
    return true;
}

//-----------------------------------------------------------------------------
void init_config()
{
    if ( ! config_input_ok() ) {
        exit(1);
    }
    int ccs,cc;
    typia=malloc(nion*sizeof*typia);
    cc=0;
    for ( int it=0; it <= ntype-1;it++) {
        invenionit[it]=1./nionit[it];
        ccs = cc;
        cc  = cc + nionit[it];
        for ( int ia=ccs;ia<=(cc-1);ia++) {
            /***********************************/
            /*        global i                 */
            /***********************************/
            typia[ia]      = it;
            massia[ia]     = massit[it];
            invemassia[ia] = 1.0/massit[it];
            /***********************************/
            /*       coulombic field           */
            /***********************************/
            if ( lcoul ) {
                qia[ia]        = qit[it]; 
                for (int i=0;i<3;i++){     
                    dipia[ia][i] = dipit[it][i];  
                    for (int j=0;j<3;j++){ 
                        quadia[ia][i][j]=quadit[it][i][j];
                    }
                }
            }
            /***********************************/
        }
    }
    nionit[ntype]=nion;
    invenionit[ntype]=1./nionit[ntype];
    info_config();

}
//-----------------------------------------------------------------------------
void info_config(){

    if (ionode) {  
        printf("configname                      : %s\n",configname);
        printf("number of types  (ntype)        : %d\n",ntype);
        printf("number of ions   (nion)         : %d\n",nion);
        putchar('\n');
    }
    info_simuCell();
}

//-----------------------------------------------------------------------------
void sample_config(int key){
    if (key==0){ /* first 10 ions */
        if (ionode ) {
            LSEPARATOR;
            printf("              Configuration sample (first 10)\n");
            LSEPARATOR;
            printf("   ia atype              rx              vx              fx\n");
            for (int ia=0;ia<4;ia++){
                printf("%5d %5s "ee3"\n",ia,atypia[ia],rx[ia],vx[ia],fx[ia]);
            }
            printf("   ia atype              ry              vy              fy\n");
            for (int ia=4;ia<8;ia++){
                printf("%5d %5s "ee3"\n",ia,atypia[ia],ry[ia],vy[ia],fy[ia]);
            }
            printf("   ia atype              rz              vz              fz\n");
            for (int ia=8;ia<12;ia++){
                printf("%5d %5s "ee3"\n",ia,atypia[ia],rz[ia],vz[ia],fz[ia]);
            }
            LSEPARATOR;
        }

    }

}

/*-----------------------------------------------------------------------------*/
/*
256
generated_by_gen_cubic.py
     20.99386800      0.00000000      0.00000000
      0.00000000     20.99386800      0.00000000
      0.00000000      0.00000000     20.99386800
1
A
256
Direct
A 0.0 0.0 0.0
...
B 0.5 0.5 0.5 
*/
int write_config(){

    double *xxx, *yyy, *zzz;
    xxx= malloc(nion*sizeof(*xxx));
    yyy= malloc(nion*sizeof(*yyy));
    zzz= malloc(nion*sizeof(*zzz));
    for(int ia=0;ia<nion;ia++){
        xxx[ia]=rx[ia];yyy[ia]=ry[ia];zzz[ia]=rz[ia];
    }

    /***************************************/ 
    /* Periodic Boundaries Conditions      */
    /***************************************/ 
    kardir ( nion , xxx , yyy , zzz , simuCell.B ) ;
    replace_pbc(nion,xxx,yyy,zzz);
    dirkar ( nion , xxx , yyy , zzz , simuCell.A ) ;

    FILE * fp;
    fp = fopen ("CONTFF", "w");
    if (NULL == fp) {
        perror("opening input file");
        return -1;
    }
    fprintf(fp,"%d\n",nion); 
    fprintf(fp,"%s\n",configname); 
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            fprintf(fp,FF,simuCell.A[i][j]);
        }
        fprintf(fp,"\n");
    }
    fprintf(fp,"%d\n",ntype);
    for(int it=0;it<ntype;it++){
        fprintf(fp,"%s",atypit[it]);
    }
    fprintf(fp,"\n");
    for(int it=0;it<ntype;it++){
        fprintf(fp,"%d",nionit[it]);
    }
    fprintf(fp,"\n");
    fprintf(fp,"Cartesian\n");
    for(int ia=0;ia<nion;ia++){
        fprintf(fp,"%s ",atypia[ia]);
        fprintf(fp,EE3,xxx[ia],yyy[ia],zzz[ia]);
        fprintf(fp,EE3,vx[ia],vy[ia],vz[ia]);
        fprintf(fp,EE3,fx[ia],fy[ia],fz[ia]);
        fprintf(fp,"\n");
    }
    fclose(fp);
    return 0;
    free(xxx);free(yyy);free(zzz);
}
/*-----------------------------------------------------------------------------*/
int read_config()
{

    pre_alloc_config();

    char cpos[MAX_LEN+1];
    FILE * fp;

    fp = fopen ("POSFF","r");

    if (NULL == fp ) {
        pError("POSFF not found\n");
#ifdef MPI
        MPI_Finalize();
#endif
        exit(-1);
    }
    // print out info to stdout
    io_node printf("reading configuration from file POSFF\n");
    
    // reading nion number of ions in POSFF 
    fscanf(fp, "%d", &nion);
    // reading config name
    fscanf(fp, "%s", configname);
    // reading cell parameters
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            fscanf(fp, "%lf",&simuCell.A[i][j]);
        }
    }
    lattice(&simuCell);
    rhoN=(double)nion * simuCell.inveOmega;

    // ------------------
    // type informations 
    // ------------------
    //ntype : number of types in POSFF
    fscanf(fp, "%d", &ntype); //ntype

    //atypit : types char 
    io_node printf("found type information on POSFF");
    for (int it=0;it<ntype;it++) {
        fscanf(fp,"%s",atypit[it]);
        io_node printf("%5s",atypit[it]);
    }
    //nionit : ions per type
    io_node printf("\n                                 ");
    for (int it=0;it<ntype;it++) {
        fscanf(fp,"%d",&nionit[it]);
        io_node printf("%5d",nionit[it]);
    }
    io_node putchar('\n');
    //can be call only when nion and ntype 
    //are known and  natmi,atypei readed
    alloc_config();

    // Direct or Cartesian
    fscanf(fp, "%s", cpos);
    //strcpy(cpos,buffer);

    // reading positions of ions 
    for (int ia=0;ia<nion;ia++) {
        fscanf(fp," %s ",atypia[ia]);
        fscanf(fp,"%lf %lf %lf ",&rx[ia],&ry[ia],&rz[ia]);
        if (Fposff>0) {
            fscanf(fp,"%lf %lf %lf ",&vx[ia],&vy[ia],&vz[ia]);
        }
        if (Fposff>1) {
            fscanf(fp,"%lf %lf %lf ",&fx[ia],&fy[ia],&fz[ia]);
        }
        //printf("ia %d atypia %s"ee9"\n",ia,atypia[ia],rx[ia],ry[ia],rz[ia],vx[ia],vy[ia],vz[ia],fx[ia],fy[ia],fz[ia]);
    }
    // if position are in direct coordinates => cartesian coordinates
    if ((strcmp(cpos,"Direct") == 0 ) || (strcmp(cpos,"D") == 0 )) {
        io_node printf("\natomic positions in Direct coordinates\n");
        dirkar(nion, rx, ry, rz , simuCell.A);
    }
    else{
        io_node printf("\natomic positions in Cartesian coordinates\n");
    }
    io_node putchar('\n');

    //closing POSFF
    if (fclose(fp))     { 
       io_node printf("error closing file."); 
#ifdef MPI
        MPI_Finalize();
#endif
       return -1; 
    }

    kardir(nion,rx,ry,rz,simuCell.B);
    replace_pbc(nion,rx,ry,rz); 
    dirkar(nion,rx,ry,rz,simuCell.A);
    return 0;
}


