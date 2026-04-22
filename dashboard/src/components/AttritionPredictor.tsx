import React, { useState } from 'react'

export const AttritionPredictor: React.FC = () => {
  const [formData, setFormData] = useState({
    Age: 30,
    Department: 'Research & Development',
    JobRole: 'Research Scientist',
    MonthlyIncome: 5000,
    JobSatisfaction: 3,
    PerformanceRating: 3,
    OverTime: 'No',
    YearsAtCompany: 5,
    WorkLifeBalance: 3,
    EngagementScore: 4.0
  })

  const [prediction, setPrediction] = useState<any>(null)
  const [loading, setLoading] = useState(false)

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault()
    setLoading(true)
    try {
      const response = await fetch('http://localhost:8000/predict-attrition', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(formData)
      })
      const data = await response.json()
      setPrediction(data)
    } catch (err) {
      console.error('Prediction failed:', err)
      alert('Failed to connect to backend. Ensure People Ops AI backend is running on port 8000.')
    } finally {
      setLoading(false)
    }
  }

  return (
    <div className="attrition-module">
      <div className="module-header">
        <h1>Retention Risk AI</h1>
        <p>Analyze employee attrition probability using machine learning fingerprints.</p>
      </div>

      <div className="g2">
        <form className="card" onSubmit={handleSubmit}>
          <h3><span className="ico">📝</span> Employee Features</h3>
          <div className="form-grid">
            <div className="f-group">
              <label>Age</label>
              <input type="number" value={formData.Age} onChange={e => setFormData({...formData, Age: parseInt(e.target.value)})} />
            </div>
            <div className="f-group">
              <label>Monthly Income ($)</label>
              <input type="number" value={formData.MonthlyIncome} onChange={e => setFormData({...formData, MonthlyIncome: parseFloat(e.target.value)})} />
            </div>
            <div className="f-group">
              <label>Department</label>
              <select value={formData.Department} onChange={e => setFormData({...formData, Department: e.target.value})}>
                <option>Sales</option>
                <option>Research & Development</option>
                <option>Human Resources</option>
              </select>
            </div>
            <div className="f-group">
              <label>Overtime</label>
              <select value={formData.OverTime} onChange={e => setFormData({...formData, OverTime: e.target.value})}>
                <option>Yes</option>
                <option>No</option>
              </select>
            </div>
            <div className="f-group">
              <label>Job Satisfaction (1-4)</label>
              <input type="range" min="1" max="4" value={formData.JobSatisfaction} onChange={e => setFormData({...formData, JobSatisfaction: parseInt(e.target.value)})} />
            </div>
            <div className="f-group">
              <label>Work-Life Balance (1-4)</label>
              <input type="range" min="1" max="4" value={formData.WorkLifeBalance} onChange={e => setFormData({...formData, WorkLifeBalance: parseInt(e.target.value)})} />
            </div>
          </div>
          <button type="submit" className="primary-btn full" disabled={loading}>
            {loading ? 'Analyzing...' : 'Predict Retention Risk'}
          </button>
        </form>

        <div className="card prediction-result">
          <h3><span className="ico">🔮</span> AI Analysis Result</h3>
          {!prediction ? (
            <div className="empty-state">
              <div className="empty-icon">🤖</div>
              <p>Waiting for data input. Submit the form to generate a risk profile.</p>
            </div>
          ) : (
            <div className="result-content">
              <div className={`risk-badge large ${prediction.risk_category.toLowerCase()}`}>
                {prediction.risk_category} RISK
              </div>
              <div className="prob-display">
                <span className="prob-val">{(prediction.attrition_probability * 100).toFixed(1)}%</span>
                <span className="prob-lbl">Probability of Attrition</span>
              </div>
              <div className="model-info">
                <span>Model: {prediction.model_used}</span>
              </div>
            </div>
          )}
        </div>
      </div>

      <style>{`
        .form-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; margin-bottom: 24px; }
        .f-group { display: flex; flex-direction: column; gap: 6px; }
        .f-group label { font-size: 0.75rem; font-weight: 600; color: var(--text-grey); }
        .f-group input, .f-group select { 
          padding: 10px; 
          border: 1px solid var(--border-light); 
          border-radius: 8px; 
          background: #f9fafb; 
          font-family: inherit; 
          font-size: 0.85rem;
        }
        .primary-btn.full { width: 100%; padding: 12px; margin-top: 8px; }

        .prediction-result { display: flex; flex-direction: column; }
        .empty-state { flex: 1; display: grid; place-content: center; text-align: center; color: var(--text-grey); }
        .empty-icon { font-size: 3rem; margin-bottom: 12px; opacity: 0.3; }
        
        .result-content { 
          flex: 1; 
          display: flex; 
          flex-direction: column; 
          align-items: center; 
          justify-content: center; 
          gap: 20px;
          animation: fadeIn 0.5s ease-out;
        }
        .risk-badge.large { 
          padding: 8px 24px; 
          border-radius: 30px; 
          font-size: 1.2rem; 
          font-weight: 800; 
        }
        .risk-badge.low { background: #d1fae5; color: #065f46; }
        .risk-badge.medium { background: #fef3c7; color: #92400e; }
        .risk-badge.high { background: #fee2e2; color: #991b1b; }
        
        .prob-display { text-align: center; }
        .prob-val { font-size: 4rem; font-weight: 900; color: var(--text-main); display: block; line-height: 1; }
        .prob-lbl { font-size: 0.9rem; color: var(--text-grey); font-weight: 500; }
        
        .model-info { font-size: 0.7rem; color: var(--text3); }

        @keyframes fadeIn { from { opacity: 0; transform: scale(0.95); } to { opacity: 1; transform: scale(1); } }
      `}</style>
    </div>
  )
}
