#ifndef JALLEGRO_MYPANEL_H
#define JALLEGRO_MYPANEL_H

#include "jallegro/Panel.h"
#include <iostream>
using namespace jallegro;

class GameOfLifePanel : public Panel {

    private:
        int* world;
        int n;
        int m;
        int p_x;
        int p_y;
        int nthx;

    public:
        #define k(r,c) (r * m + c)
        GameOfLifePanel(int* world, int n, int m, int p_x, int p_y, int nthx){
            this -> world = world;
            this -> n = n;
            this -> m = m;
            this -> p_x = p_x;
            this -> p_y = p_y;
            this -> nthx = nthx;
        }

    protected:
       void paintComponent(Graphics g){
            
            Panel::paintComponent(g);
            int blockSize = 10;
          
            for(int i = 0; i < n; i++){
                for(int j = 0; j < m; j++){
                   
                    // Impostazione del colore di riempimento e disegno del rettangolo
                    if (world[k(i, j)] == 1) {
                        g.setColor(Color::blue());
                    } else {
                        g.setColor(Color::black());
                    }
                    g.fillRect(i * blockSize, j * blockSize, blockSize, blockSize);

                    if(world[k(i, j)] == 0) {
                        // Impostazione del colore del bordo e disegno di un rettangolo più piccolo sopra di esso
                        g.setColor(Color::white()); // Imposta il colore del bordo
                        g.fillRect(i * blockSize + 1, j * blockSize + 1, blockSize - 2, blockSize - 2);
                    }
                }
            }
        }
};

#endif 