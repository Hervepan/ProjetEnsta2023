#include <SFML/Window/Keyboard.hpp>
#include <ios>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <chrono>
#include <mpi.h>
#include "cartesian_grid_of_speed.hpp"
#include "vortex.hpp"
#include "cloud_of_points.hpp"
#include "runge_kutta.hpp"
#include "screen.hpp"


using namespace std;

auto readConfigFile( std::ifstream& input )
{
    using point=Simulation::Vortices::point;

    int isMobile;
    std::size_t nbVortices;
    Numeric::CartesianGridOfSpeed cartesianGrid;
    Geometry::CloudOfPoints cloudOfPoints;
    constexpr std::size_t maxBuffer = 8192;
    char buffer[maxBuffer];
    std::string sbuffer;
    std::stringstream ibuffer;
    // Lit la première ligne de commentaire :
    input.getline(buffer, maxBuffer); // Relit un commentaire
    input.getline(buffer, maxBuffer);// Lecture de la grille cartésienne
    sbuffer = std::string(buffer,maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    double xleft, ybot, h;
    std::size_t nx, ny;
    ibuffer >> xleft >> ybot >> nx >> ny >> h;
    cartesianGrid = Numeric::CartesianGridOfSpeed({nx,ny}, point{xleft,ybot}, h);
    input.getline(buffer, maxBuffer); // Relit un commentaire
    input.getline(buffer, maxBuffer); // Lit mode de génération des particules
    sbuffer = std::string(buffer,maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    int modeGeneration;
    ibuffer >> modeGeneration;
    if (modeGeneration == 0) // Génération sur toute la grille 
    {
        std::size_t nbPoints;
        ibuffer >> nbPoints;
        cloudOfPoints = Geometry::generatePointsIn(nbPoints, {cartesianGrid.getLeftBottomVertex(), cartesianGrid.getRightTopVertex()});
    }
    else 
    {
        std::size_t nbPoints;
        double xl, xr, yb, yt;
        ibuffer >> xl >> yb >> xr >> yt >> nbPoints;
        cloudOfPoints = Geometry::generatePointsIn(nbPoints, {point{xl,yb}, point{xr,yt}});
    }
    // Lit le nombre de vortex :
    input.getline(buffer, maxBuffer); // Relit un commentaire
    input.getline(buffer, maxBuffer); // Lit le nombre de vortex
    sbuffer = std::string(buffer, maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    try {
        ibuffer >> nbVortices;        
    } catch(std::ios_base::failure& err)
    {
        std::cout << "Error " << err.what() << " found" << std::endl;
        std::cout << "Read line : " << sbuffer << std::endl;
        throw err;
    }
    Simulation::Vortices vortices(nbVortices, {cartesianGrid.getLeftBottomVertex(),
                                               cartesianGrid.getRightTopVertex()});
    input.getline(buffer, maxBuffer);// Relit un commentaire
    for (std::size_t iVortex=0; iVortex<nbVortices; ++iVortex)
    {
        input.getline(buffer, maxBuffer);
        double x,y,force;
        std::string sbuffer(buffer, maxBuffer);
        std::stringstream ibuffer(sbuffer);
        ibuffer >> x >> y >> force;
        vortices.setVortex(iVortex, point{x,y}, force);
    }
    input.getline(buffer, maxBuffer);// Relit un commentaire
    input.getline(buffer, maxBuffer);// Lit le mode de déplacement des vortex
    sbuffer = std::string(buffer,maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    ibuffer >> isMobile;
    return std::make_tuple(vortices, isMobile, cartesianGrid, cloudOfPoints);
}


int main( int nargs, char* argv[] )
{
    MPI_Init(&nargs, &argv);
    int rank, size, calcSize, calcRank;
    MPI_Comm global;
    MPI_Comm calcul;
    MPI_Comm_dup(MPI_COMM_WORLD,&global);
    MPI_Comm_rank(global, &rank);
    MPI_Comm_size(global, &size);

    int color = (rank == 0) ? 0 : 1;

    //Create a new communicator for the process that do calculation
    MPI_Comm_split(global, color, rank - 1, &calcul);
    MPI_Comm_size(calcul, &calcSize);
    MPI_Comm_rank(calcul, &calcRank);
    MPI_Request request = MPI_REQUEST_NULL;
    MPI_Status status;

    char const* filename;
    if (nargs==1)
    {
        std::cout << "Usage : vortexsimulator <nom fichier configuration>" << std::endl;
        return EXIT_FAILURE;
    }

    filename = argv[1];
    std::ifstream fich(filename);
    auto config = readConfigFile(fich);
    fich.close();

    std::size_t resx=800, resy=600;
    if (nargs>3)
    {
        resx = std::stoull(argv[2]);
        resy = std::stoull(argv[3]);
    }



    auto vortices = std::get<0>(config);
    auto isMobile = std::get<1>(config);
    auto grid     = std::get<2>(config);
    auto cloud    = std::get<3>(config);

    grid.updateVelocityField(vortices);

    int gridSize = 2*std::get<0>(grid.cellGeometry())*std::get<1>(grid.cellGeometry()); // Each cell contains the speed vector (x2)
    int vorticesSize = vortices.numberOfVortices()*3; // a vortices is defined by coordinates + intensity (x3)
    int cloudSize = cloud.numberOfPoints()*2; // a point is two coordinates (x2)
    
    std::vector<int> batchSizes(calcSize, cloudSize/calcSize);

    //We calculate the size of the uneven batch in case the number of points is not divisible by the number of process 
    int unevenBatch= cloudSize - ((int) cloudSize/calcSize)*(calcSize - 1);
    batchSizes[0] = unevenBatch;


    //We have to store the displacement because of the unevenBatch we can't use scatter 
    std::vector<int> displs(calcSize, 0);
    for (int i = 1; i < size - 1; i ++) {
        displs[i] = displs[i - 1] + batchSizes[i - 1];
    }

    grid.updateVelocityField(vortices);
    bool advance = false;
    bool animate=false;

    //Running is here to coordinates all the process and to know when we have to stop
    bool running=true;
    double dt = 0.1;

    if (rank == 0) {
        std::cout << "######## Vortex simultor ########" << std::endl << std::endl;
        std::cout << "Press P for play animation " << std::endl;
        std::cout << "Press S to stop animation" << std::endl;
        std::cout << "Press right cursor to advance step by step in time" << std::endl;
        std::cout << "Press down cursor to halve the time step" << std::endl;
        std::cout << "Press up cursor to double the time step" << std::endl;

        Graphisme::Screen myScreen( {resx,resy}, {grid.getLeftBottomVertex(), grid.getRightTopVertex()} );
        while (running) {
            auto start = std::chrono::system_clock::now();

            bool advance = false;

            // on inspecte tous les évènements de la fenêtre qui ont été émis depuis la précédente itération
            sf::Event event;
            while (myScreen.pollEvent(event)) {
                bool send = false;

                // évènement "fermeture demandée" : on ferme la fenêtre
                if (event.type == sf::Event::Closed) {
                    running = false;
                    send = true;
                }

                if (event.type == sf::Event::Resized) {
                    // on met à jour la vue, avec la nouvelle taille de la fenêtre
                    myScreen.resize(event);
                }

                if (sf::Keyboard::isKeyPressed(sf::Keyboard::P) && !animate) {
                    animate = true;
                    send = true;
                }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::S) && animate) {
                    animate = false;
                    send = true;
                }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up)) {
                    dt *= 2;
                    send = true;
                }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down)) {
                    dt /= 2;
                    send = true;
                }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right) && !advance) {
                    advance = true;
                    send = true;
                }

                //If send is true, we send all the variable needed for calculation to the process 1 
                if (send)
                {
                    MPI_Isend(&animate, 1, MPI_CXX_BOOL, 1, 0, global, &request);
                    MPI_Isend(&advance, 1, MPI_CXX_BOOL, 1, 0, global, &request);
                    MPI_Isend(&running, 1, MPI_CXX_BOOL, 1, 0, global, &request);
                    MPI_Isend(&dt, 1, MPI_DOUBLE, 1, 0, global, &request);
                }

                if (!running)
                {
                    myScreen.close();
                    break;
                }
            }
            if (animate | advance) {
                // We get back the state calculated by the other process. If the vortices aren't moving we don't need to get the grid or the vortices data 

                if (isMobile) {
                    MPI_Bcast(grid.data(), gridSize, MPI_DOUBLE, 1, global);
                    MPI_Bcast(vortices.data(), vorticesSize, MPI_DOUBLE, 1, global);
                    MPI_Recv(cloud.data(), cloudSize, MPI_DOUBLE, 1, 0, global, &status);
                } else {
                    MPI_Recv(cloud.data(), cloudSize, MPI_DOUBLE, 1, 0, global, &status);
                }
            }
            myScreen.clear(sf::Color::Black);
            std::string strDt = std::string("Time step : ") + std::to_string(dt);
            myScreen.drawText(strDt, Geometry::Point<double>{50, double(myScreen.getGeometry().second-96)});
            myScreen.displayVelocityField(grid, vortices);
            myScreen.displayParticles(grid, vortices, cloud);
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> diff = end - start;
            std::string str_fps = std::string("FPS : ") + std::to_string(1./diff.count());
            myScreen.drawText(str_fps, Geometry::Point<double>{300, double(myScreen.getGeometry().second-96)});
            myScreen.display();
        }
    }
    // He will take care of dispatching the work 
    else if (rank == 1)
    {


        while (running) {
            int flag = 0;
            //We receive the informations that we need from process 0
            MPI_Iprobe(0, 0, global, &flag, &status);
            if (flag)
            {
                MPI_Recv(&animate, 1, MPI_CXX_BOOL, 0, 0, global, &status);
                MPI_Recv(&advance, 1, MPI_CXX_BOOL, 0, 0, global, &status);
                MPI_Recv(&running, 1, MPI_CXX_BOOL, 0, 0, global, &status);
                MPI_Recv(&dt, 1, MPI_DOUBLE, 0, 0, global, &status);

                if (!running)
                    break;
            }

            if (animate | advance) {
                //We send to all the other process that are in the calcul communicator
                MPI_Bcast(&running, 1, MPI_CXX_BOOL, 0, calcul);
                MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, calcul);

                //get the size of the local batch
                int localSize = batchSizes[calcRank];

                Geometry::CloudOfPoints local_cloud(localSize);

                //Scatter everything to all the other processes 
                MPI_Scatterv(cloud.data(), batchSizes.data(), displs.data(), MPI_DOUBLE,
                             local_cloud.data(), localSize, MPI_DOUBLE, 0, calcul);

                Geometry::CloudOfPoints new_local_cloud(localSize);

                if (isMobile) {
                    new_local_cloud = Numeric::solve_RK4_movable_vortices(dt, grid, vortices, local_cloud);
                } else {
                    new_local_cloud = Numeric::solve_RK4_fixed_vortices(dt, grid, local_cloud);
                }

                MPI_Gatherv(new_local_cloud.data(), localSize, MPI_DOUBLE, cloud.data(), batchSizes.data(),
                            displs.data(), MPI_DOUBLE, 0, calcul);

                // send back grid, vortices (to everyone) and cloud to the display
                if (isMobile) {
                    MPI_Bcast(grid.data(), gridSize, MPI_DOUBLE,  1, global);
                    MPI_Bcast(vortices.data(), vorticesSize, MPI_DOUBLE, 1, global);
                    MPI_Isend(cloud.data(), cloudSize, MPI_DOUBLE, 0, 0, global, &request);
                } else {
                    MPI_Isend(cloud.data(), cloudSize, MPI_DOUBLE, 0, 0, global, &request);
                }
            }
        }
    }else{

        while (running)
        {
            MPI_Bcast(&running, 1, MPI_CXX_BOOL, 0, calcul);
            if (!running)
                break;
            MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, calcul);

            // receive particules
            int localSize = batchSizes[calcRank];
            Geometry::CloudOfPoints local_cloud(localSize);
            MPI_Scatterv(cloud.data(), batchSizes.data(), displs.data(), MPI_DOUBLE,
                         local_cloud.data(), localSize, MPI_DOUBLE, 0, calcul);

            // compute particules
            auto new_local_cloud = Numeric::solve_RK4_fixed_vortices(dt, grid, local_cloud);

            // send back particules
            MPI_Gatherv(new_local_cloud.data(), localSize, MPI_DOUBLE, cloud.data(), batchSizes.data(),
                        displs.data(), MPI_DOUBLE, 0, calcul);

            // get updated grid and vortices
            if (isMobile) {
                MPI_Bcast(grid.data(), gridSize, MPI_DOUBLE, 1, global);
                MPI_Bcast(vortices.data(), vorticesSize, MPI_DOUBLE, 1, global);
            }
        }
    }
        
    MPI_Finalize();
    return EXIT_SUCCESS;
 }