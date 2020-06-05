add_Library(processes abstract_process.h
            adsorption.h 
            adsorption.cpp
            diffusion.h
            diffusion.cpp
            factory_process.h
            factory_process.cpp
            parameters.h
            process.h
            process.cpp
            SurfaceReaction.h
            SurfaceReaction.cpp
        )

target_link_libraries (processes essentials)
