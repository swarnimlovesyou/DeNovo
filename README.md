# MedToXAi Platform

AI-powered molecular toxicity prediction system with OCR-based ingredient extraction and intelligent chemical safety assessment.

---

## Overview

MedToXAi is an end-to-end AI platform designed to analyze pharmaceutical products and predict molecular toxicity across multiple biological endpoints. The system combines OCR-based image analysis, cheminformatics, machine learning models, and large language models to support early-stage drug safety assessment.

---

## Features

- Image-based medicine label analysis using OCR
- Automated ingredient extraction from uploaded images
- Molecular toxicity prediction across 12 biological endpoints
- AI-assisted chemical safety interpretation using Groq LLaMA 3.3
- Real-time analytical visualizations
- High-performance backend with caching and rate limiting

---

## Quick Start

### Prerequisites
- Python 3.8 or higher
- Node.js 16 or higher
- Groq API key
- Supabase account (optional)

### Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/medtox-scan-ai.git
cd medtox-scan-ai

# Backend setup
cd backend
cp .env.example .env
# Configure environment variables
pip install -r requirements.txt

# Frontend setup
cd ../frontend
npm install

# Run the platform
# Backend: python app.py (port 5000)
# Frontend: npm start (port 3000)
Project Structure
pgsql
Copy code
medtox-scan-ai/
├── backend/           Flask API server
├── frontend/          React application
├── model-training/    Machine learning training pipeline
├── docs/              Documentation
└── tests/             Automated tests
Documentation
Quick Start Guide: docs/guides/QUICK_START.md

Deployment Guide: docs/guides/DEPLOYMENT_GUIDE.md

Model Training Guide: docs/training/MODEL_TRAINING_GUIDE.md

API Documentation: docs/guides/API_DOCUMENTATION.md

Full Documentation Index: docs/README.md

Model Training
Train toxicity prediction models using curated datasets:

bash
Copy code
cd model-training
python scripts/train_models.py --data data/tox21_data.csv
Refer to the Model Training Guide for configuration and evaluation details.

Technologies
Backend
Flask for API services

XGBoost for toxicity prediction models

RDKit for molecular descriptor computation

Groq AI (LLaMA 3.3) for chemical reasoning

Supabase (PostgreSQL) for data storage

Frontend
React 18 for user interface

Tailwind CSS for styling

Recharts for data visualization

Tesseract.js for OCR processing

React Hot Toast for notifications

Performance
Average prediction time: ~200 ms per molecule

Batch processing throughput: ~50 molecules per minute

Mean ROC-AUC: >84%

API response latency: 500–800 ms

Cache hit rate: approximately 60–70%

Security
Input validation and sanitization

API rate limiting (default: 60 requests per minute)

CORS protection

Secure environment variable handling

Encrypted database connections

Contributing
Contributions are welcome. Please review the Contributing Guide before submitting pull requests.

License
This project is released under the MIT License. See the LICENSE file for details.

Acknowledgments
Tox21 Challenge for toxicity datasets

RDKit for cheminformatics tooling

Groq for LLM inference

Supabase for database infrastructure

Support
Documentation: docs/README.md

Issue tracking: GitHub Issues

Contact: support@medtoxai.com
