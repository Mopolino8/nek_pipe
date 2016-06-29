#include <stdio.h>

extern int wrap_revolve(int* check,int* capo,int* fine,int *snaps_in,int* info); 

int main(int argc, char* argv[]) {
  int check, capo, oldcapo, fine, snaps_in, info, act,whattodo;
  check=-1; capo=0; fine=100;
  snaps_in=5;
  info=2;
  act=0;
  do {
    oldcapo=capo;
    whattodo=wrap_revolve(&check, &capo, &fine, &snaps_in, &info);
    switch(whattodo) {
      case 1: 
        printf("Advance from %d to %d.\n", oldcapo, capo);
        break;
      case 2:
        printf("Store in checkpoint number %d\n",check);
        break;
      case 3: 
        printf("First turn: Initialize adjoints and reverse first step.\n");
        break;
      case 4: 
        printf("Forward and reverse one step.\n");
        break;
      case 5: 
        printf("Restore in checkpoint number %d\n",check);
        break;
      case -1:
        printf("Error!");
        break;
    }
  } while(whattodo != 6);


    

}
