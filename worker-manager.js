class WorkerPool {
    constructor(workerScript, size = 0) {
        this.size = size || (navigator.hardwareConcurrency || 4);
        this.workers = [];
        this.queue = [];
        this.activeTasks = new Map(); // worker -> task
        this.workerScript = workerScript;
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
        // Cancel existing? The user asked for "Cancel in-flight computations"
        // This pool is generic, but we can clear the queue.
        // For strict cancellation, we should terminate workers or ignore results.
        // Let's implement strict cancellation: clear queue and restart workers if busy.
        
        this.cancelAll();

        const chunks = Math.ceil(totalRays / chunkSize);
        let completedChunks = 0;
        
        // Generate request ID to track this specific job
        const jobId = Date.now();

        for (let i = 0; i < chunks; i++) {
            const start = i * chunkSize;
            const count = Math.min(chunkSize, totalRays - start);
            
            const taskParams = { ...params, startIndex: start, count: count, requestId: jobId };
            
            this.queue.push({
                params: taskParams,
                onResult: (msg) => {
                    // Ignore if stale (though cancelAll should prevent this)
                    if (msg.requestId !== jobId) return;
                    
                    if (onPartialResult) onPartialResult(msg.data);
                    
                    completedChunks++;
                    if (completedChunks === chunks && onComplete) {
                        onComplete();
                    }
                }
            });
        }
        
        this.processQueue();
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
