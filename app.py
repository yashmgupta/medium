# app.py

import streamlit as st
from data import research_summaries
from collections import defaultdict

# Set page configuration
st.set_page_config(
    page_title="📚 Research Reels",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom CSS for better styling
def local_css(file_name):
    with open(file_name) as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

# Uncomment if you have a CSS file
# local_css("styles.css")

# Initialize session state variables
if 'current_index' not in st.session_state:
    st.session_state.current_index = 0

if 'likes' not in st.session_state:
    st.session_state.likes = defaultdict(bool)

if 'comments' not in st.session_state:
    st.session_state.comments = defaultdict(list)

if 'selected_categories' not in st.session_state:
    st.session_state.selected_categories = []

# Sidebar for category selection
st.sidebar.title("📖 Categories")

# Get unique categories
categories = sorted(list(set([item['category'] for item in research_summaries])))

# Multi-select for categories
selected_categories = st.sidebar.multiselect(
    "Select Categories",
    options=categories,
    default=categories,  # All selected by default
)

st.session_state.selected_categories = selected_categories

# Filter summaries based on selected categories
filtered_summaries = [
    summary for summary in research_summaries
    if summary['category'] in st.session_state.selected_categories
]

# Search bar
st.title("📚 Research Reels")
search_query = st.text_input("🔍 Search topics...", value="", placeholder="Enter keywords to search...")

# Function to filter based on search
def filter_summaries(summaries, query):
    if not query:
        return summaries
    query = query.lower()
    return [
        summary for summary in summaries
        if query in summary['title'].lower()
        or query in summary['summary'].lower()
        or any(query in tag.lower() for tag in summary['tags'])
    ]

# Apply search filter
if search_query:
    filtered_summaries = filter_summaries(filtered_summaries, search_query)

# Handle case with no results
if not filtered_summaries:
    st.warning(f"No results found for '{search_query}' in selected categories.")
    st.stop()

# Pagination: Determine the current summary
current_summary = filtered_summaries[st.session_state.current_index]

# Display the summary in an expandable section
with st.container():
    st.subheader(current_summary['title'])
    st.write(current_summary['summary'])
    st.write(f"**Category:** {current_summary['category']}")
    st.write(f"**Tags:** {', '.join(current_summary['tags'])}")

    # Action buttons
    col1, col2, col3 = st.columns(3)

    with col1:
        like_button = st.button(
            "❤️ Like",
            key=f"like_{current_summary['id']}",
            on_click=lambda: toggle_like(current_summary['id']),
        )
        if st.session_state.likes[current_summary['id']]:
            st.markdown("**You liked this!**")
    
    with col2:
        comment_input = st.text_input(
            "💬 Add a comment...",
            key=f"comment_input_{current_summary['id']}",
            placeholder="Type your comment here..."
        )
        if st.button("Submit Comment", key=f"submit_comment_{current_summary['id']}"):
            if comment_input:
                st.session_state.comments[current_summary['id']].append(comment_input)
                st.success("Comment submitted!")
    
    with col3:
        share_button = st.button("🔗 Share", key=f"share_{current_summary['id']}")
        if share_button:
            st.info("Share functionality is not implemented yet.")

    # Display comments if any
    if st.session_state.comments[current_summary['id']]:
        st.markdown("### 💬 Comments")
        for idx, comment in enumerate(st.session_state.comments[current_summary['id']], 1):
            st.write(f"{idx}. {comment}")

# Navigation buttons
st.markdown("---")
nav_col1, nav_col2, nav_col3 = st.columns([1, 2, 1])

with nav_col1:
    if st.button("⬅️ Previous", key="prev_button"):
        if st.session_state.current_index > 0:
            st.session_state.current_index -= 1
        else:
            st.warning("You're at the first summary.")

with nav_col2:
    st.write("")  # Spacer

with nav_col3:
    if st.button("Next ➡️", key="next_button"):
        if st.session_state.current_index < len(filtered_summaries) - 1:
            st.session_state.current_index += 1
        else:
            st.warning("You've reached the last summary.")

# Pagination Indicator
st.markdown(
    f"**Summary {st.session_state.current_index + 1} of {len(filtered_summaries)}**"
)

# Function to toggle like
def toggle_like(summary_id):
    st.session_state.likes[summary_id] = not st.session_state.likes[summary_id]
