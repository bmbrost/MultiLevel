# nested

Models for nested (or clustered) data structures (e.g., repeated measurements, grouped data, etc.), such that observations $i=1,\ldots,n_{j}$ are clustered within groups $j=1,\ldots,J$. 

- **normal:** Guassian data models for continuous observations
- **probit:** Probit data models (via data augmentation) for binary observations
- **poisson:** Poisson data models for counts

Naming conventions:
- **Varying-intercept models:** *.varying.intercept
- **Varying-intercept and varying-slopes models:** *.varying.coef
- **"Mixed effects" models that contain both fixed and random effects:** *.mixed
- **Models that include group-level predictors:** *.group
