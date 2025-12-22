// Physics Worker
importScripts('schwarzschild.js');

let isInitialized = false;

Module.onRuntimeInitialized = function() {
    isInitialized = true;
    postMessage({ type: 'INIT_COMPLETE' });
};

onmessage = function(e) {
    if (!isInitialized) return;
    
    const msg = e.data;
    
    if (msg.type === 'COMPUTE') {
        const { 
            startIndex, count, 
            observer_r, impact_min, impact_max, 
            num_theta, num_phi, num_impact, 
            use_spherical, num_rays_2d,
            lambda_step, lambda_max,
            requestId
        } = msg.params;
        
        try {
            // Call C++ function
            const resultVector = Module.compute_geodesic_batch(
                startIndex, count,
                observer_r, impact_min, impact_max,
                num_theta, num_phi, num_impact,
                use_spherical, num_rays_2d,
                lambda_step, lambda_max
            );
            
            // Copy data to Float32Array for transfer
            // Emscripten vectors need to be accessed carefully or copied
            // .size(), .get(i) are slow in loop.
            // There is a helper for this usually, but copying manually is safer for now.
            // Better: use direct memory view if possible, but vector is dynamic.
            
            const size = resultVector.size();
            const data = new Float32Array(size);
            for (let i = 0; i < size; i++) {
                data[i] = resultVector.get(i);
            }
            
            // Clean up C++ vector
            resultVector.delete();
            
            // Post result back (transferable)
            postMessage({
                type: 'RESULT',
                requestId: requestId,
                data: data,
                startIndex: startIndex,
                count: count
            }, [data.buffer]);
            
        } catch (err) {
            console.error("Worker error:", err);
            postMessage({ type: 'ERROR', message: err.toString(), requestId: requestId });
        }
    }
};
