class WorkerPool {
    constructor(workerScript, size = 0) {
        this.size = size || (navigator.hardwareConcurrency || 4);
        this.workers = [];
        this.queue = [];
        this.activeTasks = new Map(); // worker -> task
        this.workerScript = workerScript;
        this.cache = new Map(); // LRU cache for computed results
        this.init();
    }

    init() {
        console.log(`Initializing pool with ${this.size} workers...`);
        for (let i = 0; i < this.size; i++) {
            const worker = new Worker(this.workerScript);
            worker.id = i;
            worker.onmessage = (e) => this.handleMessage(worker, e.data);
            this.workers.push(worker);
        }
    }

    handleMessage(worker, msg) {
        if (msg.type === 'INIT_COMPLETE') {
            worker.ready = true;
            this.processQueue();
        } else if (msg.type === 'RESULT') {
            const task = this.activeTasks.get(worker);
            if (task && task.onResult) {
                task.onResult(msg);
            }
            // Worker is free
            this.activeTasks.delete(worker);
            this.processQueue();
        } else if (msg.type === 'ERROR') {
            console.error(`Worker ${worker.id} error:`, msg.message);
            this.activeTasks.delete(worker);
            this.processQueue();
        }
    }

    processQueue() {
        if (this.queue.length === 0) return;

        // Find idle ready workers
        for (const worker of this.workers) {
            if (worker.ready && !this.activeTasks.has(worker)) {
                const task = this.queue.shift();
                if (!task) break;

                this.activeTasks.set(worker, task);
                worker.postMessage({
                    type: 'COMPUTE',
                    params: task.params
                });
            }
        }
    }

    // Submit a job split into chunks
    submitJob(totalRays, chunkSize, params, onPartialResult, onComplete) {
        // Generate cache key
        const cacheKey = this.getCacheKey(params);
        if (this.cache.has(cacheKey)) {
            // Cache hit
            const cachedData = this.cache.get(cacheKey);
            // Move to end (LRU)
            this.cache.delete(cacheKey);
            this.cache.set(cacheKey, cachedData);
            
            // Immediate return
            if (onComplete) onComplete(cachedData); // Assume cachedData is full result?
            // Wait, we need to handle format. If cachedData is array of floats.
            // But onPartialResult usually takes a chunk.
            // If we have full data, we can just call onPartialResult with it once (or simulated chunks)
            // But let's assume onComplete handles the final update.
            // Actually, if we use cache, we skip partial updates.
            // To be safe, let's just call update with full data.
            // But onComplete signature is () -> void in my previous code?
            // No, onComplete usually doesn't take args, partial results accumulated elsewhere?
            // Let's check usage in shell.html (not visible here).
            // Usually onPartialResult is (data) => renderer.update(data).
            // So I should call onPartialResult(cachedData) then onComplete().
            
            if (onPartialResult) onPartialResult(cachedData);
            if (onComplete) onComplete();
            return;
        }

        // Cancel existing? The user asked for "Cancel in-flight computations"
        // This pool is generic, but we can clear the queue.
        // For strict cancellation, we should terminate workers or ignore results.
        // Let's implement strict cancellation: clear queue and restart workers if busy.
        
        this.cancelAll();

        const chunks = Math.ceil(totalRays / chunkSize);
        let completedChunks = 0;
        
        // Accumulate results for caching
        const fullResult = [];
        
        // Generate request ID to track this specific job
        const jobId = Date.now();
        
        // Track job progress for caching
        const jobResults = new Array(chunks);
        let jobCompleted = 0;

        for (let i = 0; i < chunks; i++) {
            const start = i * chunkSize;
            const count = Math.min(chunkSize, totalRays - start);
            
            const taskParams = { ...params, startIndex: start, count: count, requestId: jobId };
            
            this.queue.push({
                params: taskParams,
                onResult: (msg) => {
                    // Ignore if stale
                    if (msg.requestId !== jobId) return;
                    
                    if (onPartialResult) onPartialResult(msg.data);
                    
                    // Store chunk for cache
                    jobResults[i] = msg.data; // msg.data is array/Float32Array
                    jobCompleted++;
                    
                    if (jobCompleted === chunks) {
                        // Job done
                        // Assemble full result for cache
                        try {
                            // Calculate total length
                            let totalLen = 0;
                            for(const chunk of jobResults) totalLen += chunk.length;
                            
                            // Flatten
                            // If chunks are Float32Array, use set
                            let combined;
                            if (jobResults[0] instanceof Float32Array) {
                                combined = new Float32Array(totalLen);
                                let offset = 0;
                                for(const chunk of jobResults) {
                                    combined.set(chunk, offset);
                                    offset += chunk.length;
                                }
                            } else {
                                combined = [].concat(...jobResults);
                            }
                            
                            this.addToCache(cacheKey, combined);
                        } catch (e) {
                            console.warn("Failed to cache result:", e);
                        }
                        
                        if (onComplete) onComplete();
                    }
                }
            });
        }
        
        this.processQueue();
    } 
    
    // Helper for cache key
    getCacheKey(params) {
        // Round parameters to avoid cache misses on tiny float diffs
        const p = params;
        const r = (p.observer_r || 0).toFixed(1);
        const imin = (p.impact_min || 0).toFixed(1);
        const imax = (p.impact_max || 0).toFixed(1);
        const nphi = p.num_phi || 0;
        return `r:${r}_b:${imin}-${imax}_n:${nphi}`;
    }

    addToCache(key, data) {
        if (this.cache.size > 50) {
            // LRU eviction: delete first key (Map iterates in insertion order)
            const firstKey = this.cache.keys().next().value;
            this.cache.delete(firstKey);
        }
        this.cache.set(key, data);
    }

    cancelAll() {
        // Clear queue
        this.queue = [];
        
        // Terminate busy workers to stop computation immediately
        // This is aggressive but ensures "Main thread stays responsive" and "Prevents wasting CPU"
        const busyWorkers = [];
        for (const [worker, task] of this.activeTasks) {
            busyWorkers.push(worker);
        }
        
        busyWorkers.forEach(w => {
            w.terminate();
            // Replace with new worker
            const newWorker = new Worker(this.workerScript);
            newWorker.id = w.id;
            newWorker.onmessage = (e) => this.handleMessage(newWorker, e.data);
            
            // Replace in array
            const idx = this.workers.indexOf(w);
            if (idx !== -1) this.workers[idx] = newWorker;
            
            this.activeTasks.delete(w);
        });
    }
}
