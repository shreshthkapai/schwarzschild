class WorkerPool {
    constructor(workerScript, size = 0) {
        this.size = size || (navigator.hardwareConcurrency || 4);
        this.workers = [];
        this.queue = [];
        this.activeTasks = new Map();
        this.workerScript = workerScript;
        this.cache = new Map();
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

    submitJob(totalRays, chunkSize, params, onPartialResult, onComplete) {
        const cacheKey = this.getCacheKey(params);
        if (this.cache.has(cacheKey)) {
            const cachedData = this.cache.get(cacheKey);
            this.cache.delete(cacheKey);
            this.cache.set(cacheKey, cachedData);
            
            if (onPartialResult) onPartialResult(cachedData);
            if (onComplete) onComplete();
            return;
        }

        this.cancelAll();

        const chunks = Math.ceil(totalRays / chunkSize);
        const jobId = Date.now();
        const jobResults = new Array(chunks);
        let jobCompleted = 0;

        for (let i = 0; i < chunks; i++) {
            const start = i * chunkSize;
            const count = Math.min(chunkSize, totalRays - start);
            
            const taskParams = { ...params, startIndex: start, count: count, requestId: jobId };
            
            this.queue.push({
                params: taskParams,
                onResult: (msg) => {
                    if (msg.requestId !== jobId) return;
                    
                    if (onPartialResult) onPartialResult(msg.data);
                    
                    jobResults[i] = msg.data;
                    jobCompleted++;
                    
                    if (jobCompleted === chunks) {
                        try {
                            let totalLen = 0;
                            for(const chunk of jobResults) totalLen += chunk.length;
                            
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
    
    getCacheKey(params) {
        const p = params;
        const r = (p.observer_r || 0).toFixed(1);
        const imin = (p.impact_min || 0).toFixed(1);
        const imax = (p.impact_max || 0).toFixed(1);
        const nphi = p.num_phi || 0;
        return `r:${r}_b:${imin}-${imax}_n:${nphi}`;
    }

    addToCache(key, data) {
        if (this.cache.size > 50) {
            const firstKey = this.cache.keys().next().value;
            this.cache.delete(firstKey);
        }
        this.cache.set(key, data);
    }

    cancelAll() {
        this.queue = [];
        
        const busyWorkers = [];
        for (const [worker, task] of this.activeTasks) {
            busyWorkers.push(worker);
        }
        
        busyWorkers.forEach(w => {
            w.terminate();
            const newWorker = new Worker(this.workerScript);
            newWorker.id = w.id;
            newWorker.onmessage = (e) => this.handleMessage(newWorker, e.data);
            
            const idx = this.workers.indexOf(w);
            if (idx !== -1) this.workers[idx] = newWorker;
            
            this.activeTasks.delete(w);
        });
    }
}
