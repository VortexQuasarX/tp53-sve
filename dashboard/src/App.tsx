import { useState } from 'react'
import './index.css'
import { Login } from './components/Login'
import { Overview } from './components/Overview'
import { AttritionPredictor } from './components/AttritionPredictor'
import { InterviewIntelligence } from './components/InterviewIntelligence'
import { Leaderboard } from './components/Leaderboard'

type TabId = 'dashboard' | 'resume' | 'interviews' | 'leaderboard' | 'retention'

const TABS: { id: TabId; label: string; icon: string }[] = [
  { id: 'dashboard', label: 'Dashboard', icon: '🏠' },
  { id: 'resume', label: 'Resume Optimisation', icon: '📄' },
  { id: 'interviews', label: 'Interviews', icon: '🤝' },
  { id: 'leaderboard', label: 'Leaderboard', icon: '🏆' },
  { id: 'retention', label: 'Retention AI', icon: '🛡️' },
]

export default function App() {
  const [user, setUser] = useState<any>(null)
  const [tab, setTab] = useState<TabId>('dashboard')

  if (!user) {
    return <Login onLogin={setUser} />
  }

  return (
    <div className="portal-container">
      {/* Sidebar */}
      <div className="sidebar">
        <div className="sidebar-logo">
          <div className="logo-icon">EH</div>
          <span>EXPERT HIRE</span>
        </div>
        
        <div className="nav-list">
          {TABS.map(t => (
            <button 
              key={t.id} 
              className={`nav-item ${tab === t.id ? 'active' : ''}`} 
              onClick={() => setTab(t.id)}
            >
              <span className="nav-icon">{t.icon}</span>
              <span className="nav-label">{t.label}</span>
            </button>
          ))}
        </div>

        <div className="sidebar-footer">
          <div className="user-pill">
            <img src={user.avatar} alt="avatar" />
            <div className="user-info">
              <span className="user-name">{user.name}</span>
              <span className="user-role">Student</span>
            </div>
          </div>
        </div>
      </div>

      {/* Main Content */}
      <div className="main-viewport">
        <header className="main-header">
          <div className="breadcrumb">Dashboard / {TABS.find(t => t.id === tab)?.label}</div>
          <button className="primary-btn">Take Interview</button>
        </header>

        <div className="content-area">
          {tab === 'dashboard' && <Overview />}
          {tab === 'retention' && <AttritionPredictor />}
          {tab === 'interviews' && <InterviewIntelligence />}
          {tab === 'leaderboard' && <Leaderboard />}
          {/* Other tabs placeholder */}
          {['resume'].includes(tab) && (
            <div className="placeholder-view">
              <h2>{TABS.find(t => t.id === tab)?.label}</h2>
              <p>Module integration in progress...</p>
            </div>
          )}
        </div>
      </div>
    </div>
  )
}
