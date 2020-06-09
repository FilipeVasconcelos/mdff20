void print_tensorRK2(TENSOR_RK2 tens){
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            printf(ee,tens.ab[i][j]);
        }
        putchat('\n');
    }
    putchat('\n');
}
